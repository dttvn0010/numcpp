#include "index.h"
#include "np.h"
#include <cassert>
#include <memory>

Slice::Slice()
{

}

Slice::Slice(Integer start, Integer stop, Integer step)
{
    this->start = start;
    this->stop = stop;
    this->step = step;
}

Slice Slice::all()
{
    return Slice();
}

Slice Slice::head(int stop)
{
    return Slice(Integer::none(), stop, Integer::none());
}

Slice Slice::head(int stop, int step)
{
    return Slice(Integer::none(), stop, step);
}

Slice Slice::tail(int start)
{
    return Slice(start, Integer::none(), Integer::none());
}

Slice Slice::tail(int start, int step)
{
    return Slice(start, Integer::none(), step);
}

Slice Slice::step_by(int step)
{
    return Slice(Integer::none(), Integer::none(), step);
}

void Slice::resolve_index(int size) 
{
    int step = this->step.is_null ? 1 : this->step.value;
    int start = this->start.is_null? (step > 0? 0 : size-1) : this->start.value;
    int stop = this->stop.is_null? (step > 0? size : -1) : this->stop.value;

    if(start < 0) start += size;
    if(stop < 0 && !this->stop.is_null) stop += size;
    if(start < 0) start = 0;
    if(stop > size) stop = size;
    
    if((step > 0 && stop < start) || (step < 0 && stop > start)){
        stop = start;
    }

    if(step > 0){
        stop = start + ((stop - start - 1) / step) * step + step;
    }else if(step < 0){
        stop = start - ((start - stop - 1) / (-step)) * (-step) + step;
    }
    this->start = start;
    this->stop = stop;
    this->step = step;
}

int Slice::get_size() const
{
    return (stop.value - start.value) / step.value;
}
// ==================================================================================================================================
Index::Index()
{
    type = Index::INDEX_NONE;
}

void Index::copy_from(const Index& idx)
{
    type = idx.type;
    if(type == INDEX_SCALAR)
    {
        scalar = idx.scalar;
    }
    else if(type == INDEX_SLICE)
    {
        slice = std::make_unique<Slice>(*idx.slice);
    }
    else if(type == INDEX_LIST)
    {
        int *p = new int[idx.index_length];
        for(int i = 0; i < idx.index_length; i++) p[i] = idx.indexes[i];
        indexes = std::unique_ptr<int[]>(p);
        index_length = idx.index_length;
    }
}

Index::Index(const Index& idx)
{
    copy_from(idx);
}

Index& Index::operator= (const Index& idx)
{
    copy_from(idx);
    return *this;
}

Index::Index(int index)
{
    type = Index::INDEX_SCALAR;
    this->scalar = index;
}

Index::Index(const int* indexes, int n)
{
    type = Index::INDEX_LIST;
    int* p = new int[n];
    for(int i = 0; i < n; i++) p[i] = indexes[i];
    this->indexes = std::unique_ptr<int[]>(p);
    this->index_length = n;
}

Index::Index(const List<int>& indexes) : Index(indexes.get_data(), indexes.size())
{
    
}

Index::Index(const NDArray<int>& indexes) : Index(indexes.get_data(), indexes.get_shape(0))
{
    assert(indexes.get_dim() == 1);
}

Index::Index(NDArray<int>&& indexes)
{
    assert(indexes.get_dim() == 1);
    type = Index::INDEX_LIST;
    this->indexes = std::unique_ptr<int[]>(indexes.into_raw());
}

Index::Index(const boolean* indexes, int n)
{
    type = Index::INDEX_LIST;
    max_idx = n;
    index_length = 0;
    for(int i = 0; i < n; i++) if(indexes[i]) index_length += 1;
    int* p = new int[index_length];
    int idx = 0;
    for(int i = 0; i < n; i++)
    {
        if(indexes[i])
        {
            p[idx] = i;
            idx += 1;
        }
    }
    this->indexes = std::unique_ptr<int[]>(p);
}

Index::Index(const List<boolean>& indexes) : Index(indexes.get_data(), indexes.size())
{
    
}

Index::Index(const NDArray<boolean>& indexes) : Index(indexes.get_data(), indexes.get_shape(0))
{
    assert(indexes.get_dim() == 1);
}

Index::Index(const Slice& slice)
{
    type = Index::INDEX_SLICE;
    this->slice = std::make_unique<Slice>(slice);
}


bool Index::is_none() const
{
    return type == INDEX_NONE;
}

bool Index::is_scalar() const
{
    return type == INDEX_SCALAR;
}

bool Index::is_slice() const
{
    return type == INDEX_SLICE;
}

bool Index::is_list() const 
{
    return type == INDEX_LIST;
}

int Index::get_scalar() const
{
    return scalar;
}

Index Index::none()
{
    return Index();
}

Slice* Index::get_slice()
{
    return slice.get();
}

const Slice* Index::get_slice() const
{
    return slice.get();
}

int* Index::get_indexes() const
{
    return indexes.get();
}

int Index::get_max_idx() const
{
    return max_idx;
}

int Index::get_index_length() const 
{
    return index_length;
}