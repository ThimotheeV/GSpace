#pragma once // avoid redefinitions (similar as #ifdef...)
#include <memory>
#include <vector>
#include <functional>
#include <map>
#include <iostream>

int const MEM_FACTOR = 10000;

namespace unit_test
{
    struct object_memory_manager_test;
}
//Combines seg's p'ointers with indexes, allows translate from one to the other and making indexes available when their associated ptr is destroyed
template <class T>
class obj_mem_manag_c
{
    //available_objs_indexes indicate container's indexes who were free to use. available_objs contains uniqu_ptr to your object.
    std::map<T *, std::size_t> hash_ptr_to_index;
    std::vector<std::size_t> available_objs_indexes;
    std::vector<T *> in_use_objs;

    // for unit testing
    friend struct unit_test::object_memory_manager_test;

    //Custom ptr (object_ptr) with the delete fonction above
    using object_ptr = T *;

public:
    obj_mem_manag_c() = default;

    obj_mem_manag_c(const obj_mem_manag_c &) = delete;
    obj_mem_manag_c(obj_mem_manag_c &&) = delete;
    obj_mem_manag_c &operator=(const obj_mem_manag_c) = delete;
    obj_mem_manag_c &operator=(obj_mem_manag_c &&) = delete;

    auto get_new_obj()
    {
        if (available_objs_indexes.empty())
        {
            auto old_size = in_use_objs.size();
            //+20000 for the first turn
            auto new_size = (old_size + MEM_FACTOR) * 2;
            available_objs_indexes.reserve(new_size);
            in_use_objs.resize(new_size);
            //Create new index for futur uses
            for (std::size_t index = new_size - 1; index > old_size; --index)
            {
                available_objs_indexes.push_back(index);
            }
            //For the first turn if index go under 0 he go to max value of size_t
            available_objs_indexes.push_back(old_size);
        }
        auto index = available_objs_indexes.back();
        available_objs_indexes.pop_back();

        auto &ptr = in_use_objs[index];
        if (ptr == nullptr)
        {
            ptr = new T;
            hash_ptr_to_index.emplace(ptr, index);
        }

        return ptr;
    }

    //Can't call a non existante ptr
    std::size_t get_index_from_obj(object_ptr ptr) noexcept
    {
        return hash_ptr_to_index.at(ptr);
    }

    auto get_obj_from_index(std::size_t index) noexcept
    {
        return in_use_objs[index];
    }

    void free_obj(object_ptr &ptr)
    {
        auto index = get_index_from_obj(ptr);
        available_objs_indexes.push_back(index);
        (in_use_objs[index])->clean();

        ptr = nullptr;
    }

    ~obj_mem_manag_c()
    {
        for (auto ptr : in_use_objs)
        {
            delete ptr;
        }
    }
};