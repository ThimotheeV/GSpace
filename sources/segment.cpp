#include "segment.hpp"

void seg::clean()
{
    L_left_brkpt = 0;
    R_right_brkpt = 0;
    U_ident_node = 0;
    Previous_seg = nullptr;
    Next_seg = nullptr;
    Cell_lineage = nullptr;
}

indiv_c::indiv_c()
{
    First_seg_chr.resize(1);
}

indiv_c::indiv_c(int chr_nbr, node_lattice_c *node_lat) : Node_lat(node_lat)
{
    First_seg_chr = std::vector<std::array<seg *, 2>>(chr_nbr, std::array<seg *, 2>{nullptr, nullptr});
}

void indiv_c::update_homolog_seg(std::size_t chr_index, seg *homolog_seg)
{
    if (First_seg_chr[chr_index].at(1) != nullptr)
    {
        (First_seg_chr[chr_index].at(1))->set_indiv(nullptr);
    }

    if (homolog_seg != nullptr)
    {
        homolog_seg->set_indiv(this);
    }

    First_seg_chr[chr_index].at(1) = homolog_seg;
}

void indiv_c::update_indiv(std::size_t chr_index, seg *first_seg, seg *second_seg)
{
    //order because indiv_c need to have a effective first_seg
    if (first_seg == nullptr)
    {
        std::swap(first_seg, second_seg);
    }

    if (First_seg_chr[chr_index].at(0) != nullptr)
    {
        (First_seg_chr[chr_index].at(0))->set_indiv(nullptr);
    }

    if (First_seg_chr[chr_index].at(1) != nullptr)
    {
        (First_seg_chr[chr_index].at(1))->set_indiv(nullptr);
    }

    if (first_seg != nullptr)
    {
        (first_seg)->set_indiv(this);
    }

    First_seg_chr[chr_index].at(0) = first_seg;

    if (second_seg != nullptr)
    {
        (second_seg)->set_indiv(this);
    }
    First_seg_chr[chr_index].at(1) = second_seg;
}

void indiv_c::untie_chr(std::size_t chr_index, std::size_t chr)
{
    if (First_seg_chr[chr_index].at(chr) != nullptr)
    {
        (First_seg_chr[chr_index].at(chr))->set_indiv(nullptr);
        (First_seg_chr[chr_index].at(chr)) = nullptr;
    }
}

bool indiv_c::empty() const
{
    for (auto const seg : First_seg_chr)
    {
        if (seg.at(0) != nullptr)
        {
            return false;
        }
    }
    return true;
}

//indiv_stock_c//

seg *find_first_seg(seg *middle_seg)
{
    while (middle_seg->Previous_seg != nullptr)
    {
        middle_seg = middle_seg->Previous_seg;
    }
    return middle_seg;
}

std::array<seg *, 4> &random_coa_process_btw_seg(std::array<seg *, 4> &radomize_coa_order, rand_gen_c &rand_gen)
{
    //{A1, A2, B1, B2} -> {A1, B1, A2, B2}
    std::iter_swap(radomize_coa_order.begin() + 1, radomize_coa_order.begin() + 2);

    //{A1, B1, A2, B2} -> {A1, B2, A2, B1}
    if (rand_gen.rand_bool())
    {
        std::iter_swap(radomize_coa_order.begin() + 1, radomize_coa_order.begin() + 3);
    }

    return radomize_coa_order;
}

indiv_c *indiv_stock_c::new_indiv(node_lattice_c *node_lat)
{
    indiv_c *indiv;

    if (No_use.empty())
    {
        indiv = new indiv_c(Chr_nbr, node_lat);
    }
    else
    {
        indiv = No_use.back();
        indiv->Node_lat = node_lat;
        No_use.pop_back();
    }

    add_indiv(indiv);
    return indiv;
}

void indiv_stock_c::add_indiv(indiv_c *indiv)
{
    if (Stock.capacity() < (Stock.size() + 1))
    {
        Stock.reserve((Stock.size() + 1) * 2);
    }
    Stock.push_back(indiv);
    indiv->Ident = size() - 1;
}

void indiv_stock_c::clean_indiv_at_chr(std::size_t ind, std::size_t chr_index)
{
    auto indiv = Stock[ind];

    indiv->untie_chr(chr_index, 1);
    indiv->untie_chr(chr_index, 0);
}

void indiv_stock_c::erase(std::size_t ind)
{
    if (ind >= Stock.size() || ind < 0)
    {
        throw std::logic_error("( In indiv_stock_c::erase() : Try to erase a cell_node who haven't a coherant ident. Contact the developpers. I exit. )");
    }

    //Cancel overheap of erase a value in a middle of a vector
    std::iter_swap(Stock.begin() + ind, Stock.end() - 1);
    Stock[ind]->Ident = ind;

    int last_ind = Stock.size() - 1;

    for (int chr_index = 0; chr_index < Chr_nbr; ++chr_index)
    {
        clean_indiv_at_chr(last_ind, chr_index);
    }
    Stock.back()->Node_lat = nullptr;
    Stock.back()->Ident = -1;
    No_use.push_back(Stock.back());
    Stock.pop_back();
}

indiv_stock_c::~indiv_stock_c()
{
    for (auto ptr : Stock)
    {
        delete ptr;
    }
    for (auto ptr : No_use)
    {
        delete ptr;
    }
}

void indiv_stock_c::ini(lattice_c &lat, samp_param_c const &samp_param)
{
    std::size_t indiv_nbr = samp_param.n_total_sample_size;

    if (samp_param.Sample_coord_vec.size() != indiv_nbr)
    {
        throw std::logic_error("( In indiv_stock_c::ini() : All sample don't have a coord in vec. Contact the developpers. I exit. )");
    }

    auto coord_vec_itr = samp_param.Sample_coord_vec.begin();

    Stock.reserve(indiv_nbr * 10);

    for (std::size_t i = 0; i < indiv_nbr; ++i)
    {
        //Lat/node/cell intialisation
        auto node = lat.node(*coord_vec_itr++);
        //WARNING : Create a new indiv
        auto indiv_ptr = new indiv_c(Chr_nbr, node);
        add_indiv(indiv_ptr);
        add_indiv_at_node(indiv_ptr, node);
    }

    No_use.reserve(indiv_nbr * 10);
}
