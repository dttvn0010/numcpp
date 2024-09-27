#ifndef __INDEX_H__
#define __INDEX_H__
#include "base.h"
template<typename  T>
class NDArray;

struct Slice 
{
    Integer start = Integer::none();
    Integer stop = Integer::none();
    Integer step = Integer::none();

    Slice();
    Slice(Integer start, Integer stop, Integer step=Integer::none());

    static Slice all();

    static Slice head(int stop);

    static Slice head(int stop, int step);

    static Slice tail(int start);

    static Slice tail(int start, int step);

    static Slice step_by(int step);

    void resolve_index(int size);

    int get_size() const;
};

class Index
{
    const static int INDEX_NONE = 0;
    const static int INDEX_SCALAR = 1;
    const static int INDEX_SLICE = 2;
    const static int INDEX_LIST = 3;
    const static int INDEX_MASK = 4;

    int type = 0;
    int scalar = 0;
    int max_idx = -1;
    int index_length = 0;
    std::unique_ptr<Slice> slice;
    std::unique_ptr<int[]> indexes;

    void copy_from(const Index& index);
public:
    Index();

    Index(Index&&) = default;

    Index(const Index&);

    Index(int index);

    Index(const Slice& slice);

    Index(const int* indexes, int n);

    Index(const List<int>& indexes);

    Index(const NDArray<int>& indexes);

    Index(NDArray<int>&& indexes);

    Index(const boolean* indexes, int n);

    Index(const List<boolean>& indexes);

    Index(const NDArray<boolean>& indexes);

    Index& operator = (const Index&);

    Index& operator =(Index&&) = default;

    bool is_none() const;

    bool is_scalar() const;

    bool is_slice() const;

    bool is_list() const;

    int get_scalar() const;

    const Slice* get_slice() const;

    Slice* get_slice();

    int* get_indexes() const;
    
    int get_index_length() const;

    int get_max_idx() const;

    static Index none();
};
#endif