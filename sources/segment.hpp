#pragma once
#include <vector>
#include <memory>
#include <map>

#include "object_memory_manager.hpp"
#include "random_generator.hpp"
#include "migration.hpp"
#include "settings.hpp"

struct node_lattice_c;
struct indiv_c;
struct lattice_c;
struct samp_param_c;

struct seg
{
public:
  //Constructor
  seg() = default;
  //constructor by copy directly using the copy constructor of each object
  // instead of classical copy constructor that first creates the object using the empty constructor, and then copy the objects
  seg(int left_brkpt_l, int right_brkpt_r, int node_u, seg *previous_seg = nullptr, seg *next_seg = nullptr, indiv_c *cell_lineage = nullptr)
      : L_left_brkpt(left_brkpt_l),
        R_right_brkpt(right_brkpt_r),
        U_ident_node(node_u),
        Previous_seg(previous_seg),
        Next_seg(next_seg),
        Cell_lineage(cell_lineage){}; // no instructions in the function, only pre-function instructions

  //Members Data
  int L_left_brkpt{0};
  int R_right_brkpt{0};
  int U_ident_node{0};
  seg *Previous_seg{nullptr};
  seg *Next_seg{nullptr};

  indiv_c *get_indiv() const
  {
    return Cell_lineage;
  }

  void set_indiv(indiv_c *indiv)
  {
    Cell_lineage = indiv;
  }

  void clean();

private:
  //Only indiv_c can change it
  indiv_c *Cell_lineage{nullptr};
};

struct indiv_c
{
  indiv_c();
  indiv_c(int chr_nbr, node_lattice_c *node_lat);

  seg *get_chr(std::size_t chr_index, std::size_t place) const
  {
    return First_seg_chr[chr_index].at(place);
  }

  std::size_t nbr_chr() const
  {
    return First_seg_chr.size();
  }

  void update_homolog_seg(std::size_t chr_index, seg *homolog_seg);
  void update_indiv(std::size_t chr_index, seg *first_seg, seg *second_seg);
  void untie_chr(std::size_t chr_index, std::size_t chr);

  bool empty() const;

private:
  std::vector<std::array<seg *, 2>> First_seg_chr;

public:
  node_lattice_c *Node_lat{nullptr};
  int Ident{-1};
};

seg *find_first_seg(seg *middle_seg);
std::array<seg *, 4> &random_coa_process_btw_seg(std::array<seg *, 4> &radomize_coa_order, rand_gen_c &rand_gen);

struct indiv_stock_c
{
  indiv_stock_c() = default;

  indiv_stock_c(const indiv_stock_c &) = delete;
  indiv_stock_c(indiv_stock_c &&) = delete;
  indiv_stock_c &operator=(const indiv_stock_c) = delete;
  indiv_stock_c &operator=(indiv_stock_c &&) = delete;

  explicit indiv_stock_c(int chr_nbr) : Chr_nbr(chr_nbr){};
  //operator
  indiv_c *operator[](std::size_t i)
  {
    return Stock[i];
  }

  indiv_c *operator[](std::size_t i) const
  {
    return Stock[i];
  }

  auto begin()
  {
    return Stock.begin();
  }

  auto end()
  {
    return Stock.end();
  }

  auto begin() const
  {
    return Stock.cbegin();
  }

  auto end() const
  {
    return Stock.cend();
  }

  auto back()
  {
    return Stock.back();
  }
  //information
  std::size_t size()
  {
    return Stock.size();
  }

  void reserve(std::size_t i)
  {
    Stock.reserve(i);
    No_use.reserve(i);
  }

  std::size_t capacity()
  {
    return Stock.capacity();
  }

  //transformation
  indiv_c *new_indiv(node_lattice_c *node_lat);
  void add_indiv(indiv_c *indiv);
  void clean_indiv_at_chr(std::size_t ind, std::size_t chr_index);
  void erase(std::size_t ind);

  void ini(lattice_c &lat, samp_param_c const &samp_param);

  //destructor
  ~indiv_stock_c();

private:
  std::vector<indiv_c *> Stock;
  std::vector<indiv_c *> No_use;
  int Chr_nbr{1};
};
