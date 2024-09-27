#ifndef __ND_ARRAY_H__
#define __ND_ARRAY_H__
#include "base.h"
#include "index.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <math.h>

#define MAX_DIM 16

constexpr size_t NUM_OP = 17;
constexpr size_t OP_ASSIGN = 0; 
constexpr size_t OP_ADD = 1;
constexpr size_t OP_SUB = 2;
constexpr size_t OP_MUL = 3;
constexpr size_t OP_DIV = 4;
constexpr size_t OP_MOD = 5;
constexpr size_t OP_POW = 6;
constexpr size_t OP_AND = 7;
constexpr size_t OP_OR  = 8;
constexpr size_t OP_XOR = 9;
constexpr size_t OP_GT  = 10;
constexpr size_t OP_GTE = 11;
constexpr size_t OP_LT  = 12;
constexpr size_t OP_LTE = 13;
constexpr size_t OP_EQ  = 14;
constexpr size_t OP_NE  = 15;
constexpr size_t OP_MINIMUM = 16;
constexpr size_t OP_MAXIMUM = 16;

constexpr size_t NUM_REDUCE = 12;
constexpr size_t REDUCE_MAX = 1;
constexpr size_t REDUCE_MIN = 2;
constexpr size_t REDUCE_SUM = 3;
constexpr size_t REDUCE_AND = 4;
constexpr size_t REDUCE_OR  = 5;
constexpr size_t REDUCE_XOR = 6;
constexpr size_t REDUCE_MEAN = 7;
constexpr size_t REDUCE_MEDIAN = 8;
constexpr size_t REDUCE_STD = 9;
constexpr size_t REDUCE_ARGMIN = 10;
constexpr size_t REDUCE_ARGMAX = 11;
constexpr size_t REDUCE_PERCENTILE = 12;

constexpr size_t NUM_FUNC = 22;
constexpr size_t FUNC_ABS = 1;
constexpr size_t FUNC_ROUND = 2;
constexpr size_t FUNC_CEIL = 3;
constexpr size_t FUNC_FLOOR = 4;
constexpr size_t FUNC_SQUARE = 5;
constexpr size_t FUNC_SQRT = 6;
constexpr size_t FUNC_EXP = 7;
constexpr size_t FUNC_LOG = 8;
constexpr size_t FUNC_SIN = 9;
constexpr size_t FUNC_COS = 10;
constexpr size_t FUNC_TAN = 11;
constexpr size_t FUNC_ASIN = 12;
constexpr size_t FUNC_ACOS = 13;
constexpr size_t FUNC_ATAN = 14;
constexpr size_t FUNC_SINH = 15;
constexpr size_t FUNC_COSH = 16;
constexpr size_t FUNC_TANH = 17;
constexpr size_t FUNC_ASINH = 18;
constexpr size_t FUNC_ACOSH = 19;
constexpr size_t FUNC_ATANH = 20;
constexpr size_t FUNC_NEG = 21;
constexpr size_t FUNC_NOT = 22;

template<int func, typename T>
inline static T _unary_ops(const T& val)
{
    if constexpr(func == FUNC_ABS) return ::fabs(val);
    if constexpr(func == FUNC_ROUND) return ::round(val);
    if constexpr(func == FUNC_CEIL) return ::ceil(val);
    if constexpr(func == FUNC_FLOOR) return ::floor(val);
    if constexpr(func == FUNC_SQUARE) return val*val;
    if constexpr(func == FUNC_SQRT) return ::sqrt(val);
    if constexpr(func == FUNC_EXP) return ::exp(val);
    if constexpr(func == FUNC_LOG) return ::log(val);
    if constexpr(func == FUNC_SIN) return ::sin(val);
    if constexpr(func == FUNC_COS) return ::cos(val);
    if constexpr(func == FUNC_TAN) return ::tan(val);
    if constexpr(func == FUNC_ASIN) return ::asin(val);
    if constexpr(func == FUNC_ACOS) return ::acos(val);
    if constexpr(func == FUNC_ATAN) return ::atan(val);
    if constexpr(func == FUNC_SINH) return ::sinh(val);
    if constexpr(func == FUNC_COSH) return ::cosh(val);
    if constexpr(func == FUNC_TANH) return ::tanh(val);
    if constexpr(func == FUNC_ASINH) return ::asinh(val);
    if constexpr(func == FUNC_ACOSH) return ::acosh(val);
    if constexpr(func == FUNC_ATANH) return ::atanh(val);
    if constexpr(func == FUNC_NEG) return -val;
    if constexpr(func == FUNC_NOT) return ~val;
}

template<int ops, typename T, typename V>
inline static auto _bin_ops(const T& l, const V& r)
{
    if constexpr(ops == OP_ADD) return l + r;
    if constexpr(ops == OP_SUB) return l - r;
    if constexpr(ops == OP_MUL) return l * r;
    if constexpr(ops == OP_DIV) return l / r;
    if constexpr(ops == OP_MOD) return l % r;
    if constexpr(ops == OP_POW) return ::pow(l, r);
    if constexpr(ops == OP_AND) return l & r;
    if constexpr(ops == OP_OR)  return l | r;
    if constexpr(ops == OP_XOR) return l ^ r;
    if constexpr(ops == OP_GT)  return l > r;
    if constexpr(ops == OP_GTE) return l >= r;
    if constexpr(ops == OP_LT)  return l < r;
    if constexpr(ops == OP_LTE) return l <= r;
    if constexpr(ops == OP_EQ)  return l == r;
    if constexpr(ops == OP_NE)  return l != r;
    if constexpr(ops == OP_MINIMUM)  return (l < r) ? l : r;
    if constexpr(ops == OP_MAXIMUM)  return (l > r) ? l : r;
}

template<int ops, typename T, typename V>
inline static void _bin_ops_assign(T& l, const V& r) 
{
    if constexpr(ops == OP_ASSIGN) l = (T) r;
    if constexpr(ops == OP_ADD) l += r;
    if constexpr(ops == OP_SUB) l -= r;
    if constexpr(ops == OP_MUL) l *= r;
    if constexpr(ops == OP_DIV) l /= r;
    if constexpr(ops == OP_MOD) l %= r;
    if constexpr(ops == OP_AND) l &= r;
    if constexpr(ops == OP_OR)  l |= r;
    if constexpr(ops == OP_XOR) l ^= r;

}

template<int ops, typename T, typename V>
inline static auto _reduce_ops(const V& a, const T& b) 
{
    if constexpr(ops == REDUCE_SUM) return (a + b);
    if constexpr(ops == REDUCE_MAX) return (a > b ? a : b);
    if constexpr(ops == REDUCE_MIN) return (a < b ? a : b);
    if constexpr(ops == REDUCE_AND) return a & b;
    if constexpr(ops == REDUCE_OR) return a | b;
    if constexpr(ops == REDUCE_XOR) return a ^ b;
    if constexpr(ops == REDUCE_ARGMIN) return (a < b) ? 0 : 1;
    if constexpr(ops == REDUCE_ARGMAX) return (a > b) ? 0 : 1;
    if constexpr(ops == REDUCE_STD || ops == REDUCE_MEAN || ops == REDUCE_MEDIAN) return 0.0;
}

template<typename T>
void swap(T* a, T* b)
{
    T t = *a;
    *a = *b;
    *b = t;
}

template<typename T>
void swap_with_idx(T* a, T* b, int* ia, int* ib)
{
    T t = *a;
    *a = *b;
    *b = t;
    int i = *ia;
    *ia = *ib;
    *ib = i;
}

template<typename T>
int partition(T arr[], int step, int indexes[], int low, int high)
{
    // Choosing the pivot
    T pivot = arr[high*step];
 
    // Index of smaller element and indicates
    // the right position of pivot found so far
    int i = (low - 1);
 
    for (int j = low; j <= high - 1; j++) {
 
        // If current element is smaller than the pivot
        if (arr[j*step] < pivot) {
 
            // Increment index of smaller element
            i++;
            if(indexes)
            {
                swap_with_idx(&arr[i*step], &arr[j*step], &indexes[i], &indexes[j]);
            }else
            {
                swap(&arr[i*step], &arr[j*step]);
            }
        }
    }
    if(indexes)
    {
        swap_with_idx(&arr[(i+1)*step], &arr[high*step], &indexes[i+1], &indexes[high]);
    }else
    {
        swap(&arr[(i+1)*step], &arr[high*step]);
    }
    return (i + 1);
}

template<typename T>
void quick_sort(T arr[], int step, int indexes[], int low, int high)
{
    if (low < high) {
 
        // pi is partitioning index, arr[p]
        // is now at right place
        int pi = partition(arr, step, indexes, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quick_sort(arr, step, indexes, low, pi - 1);
        quick_sort(arr, step, indexes, pi + 1, high);
    }
}

template<typename T>
class Array;

template<typename T>
class Array2;

template<typename T>
class Array3;

template<typename T>
class Array4;

template<typename T>
class NDArray : public Object  {
protected:    
    const static int MIN_CONTINUOUS_CHUNK = 32;

    template<typename U> friend class NDArray;
    template<typename U> friend class Array;
    template<typename U> friend class Array2;
    template<typename U> friend class Array3;
    template<typename U> friend class Array4;

    T* _data = NULL;
    int _shapes[MAX_DIM];
    int _steps[MAX_DIM];
    int _dim;
    bool _is_view = false;
    const size_t* _p_signature = NULL;

    NDArray(T* data, int dim, const int* shapes, const int* steps, bool is_view=false, const size_t* p_signature=NULL)
    {
        assert(dim >= 0 && dim <= MAX_DIM);
        _data = data;
        _dim = dim;
        for(int i = 0; i < dim; i++)
        {
            _shapes[i] = shapes[i];
            _steps[i] = steps[i];
        }
        _is_view = is_view;
        if(is_view)
        {
            assert(p_signature != NULL);
            _signature = *p_signature;
            _p_signature = p_signature;
        }
    }

    int get_min_ax(int bound) const 
    {
        int min_ax = 0;
        int min_shape = _shapes[0];
        for(int i = 1; i <= bound; i++)
        {
            if(_shapes[i] < min_shape)
            {
                min_shape = _shapes[i];
                min_ax = i;
            }
        }
        return min_ax;
    }

    int get_best_ax() const
    {
        int bound = get_dim() - 1;
        if(bound > 1 && _shapes[bound] >= MIN_CONTINUOUS_CHUNK) bound -= 1;
        return get_min_ax(bound);
    }

    int get_best_ax(const Index* indexes, int index_size) const
    {
        int last_idx = _dim - 1;
        int step = 1;

        while(last_idx >= 0)
        {
            if(_steps[last_idx] != step) break;
            if(last_idx < index_size && 
                !(
                    indexes[last_idx].is_slice() &&
                    indexes[last_idx].get_slice()->start.value == 0 &&
                    indexes[last_idx].get_slice()->stop.value == _shapes[last_idx] &&
                    indexes[last_idx].get_slice()->step.value == 1
                )
            )
            {
                break;
            }
            step *= _shapes[last_idx];
            last_idx -= 1;
        }

        if(last_idx + 1 < _dim && _steps[last_idx + 1] < MIN_CONTINUOUS_CHUNK)
        {
            last_idx = index_size - 1;
        }
        
        int min_ax = 0;
        int min_sz = -1;

        for(int i = 0; i <= last_idx; i++)
        {
            int sz = 0;
            if(indexes[i].is_slice())
            {
                sz = indexes[i].get_slice()->get_size();
            }
            else 
            {
                sz = indexes[i].get_index_length();
            }
            if(min_sz < 0 || sz < min_sz)
            {
                min_sz = sz;
                min_ax = i;
            }
        }
        return min_ax;
    }

    template<typename U>
    static int get_best_ax(const NDArray& lhs, const NDArray<U>& rhs)
    {
        assert(lhs.get_dim() == rhs.get_dim());
        int dim = lhs.get_dim();

        int min_ax = 0;
        int min_shape = std::max(lhs._shapes[0], rhs._shapes[0]);
        for(int i = 1; i < dim; i++)
        {
            int shape_i = std::max(lhs._shapes[i], rhs._shapes[i]);
            if(shape_i < min_shape)
            {
                min_shape = shape_i;
                min_ax = i;
            }
        }
        return min_ax;
    }

    NDArray _get_slice(int ax, int idx) const
    {
        int sub_shapes[MAX_DIM];
        int sub_steps[MAX_DIM];
        int sub_dim = 0;
        for(int i = 0; i < _dim; i++)
        {
            if(i != ax)
            {
                sub_shapes[sub_dim] = _shapes[i];
                sub_steps[sub_dim] = _steps[i];
                sub_dim += 1;
            }
        }
        return NDArray(_data + idx * _steps[ax], sub_dim, sub_shapes, sub_steps, true, get_p_signature());
    }

    NDArray _get_slice(int ax, int start, int end, int step=1) const
    {
        int sub_shapes[MAX_DIM];
        int sub_steps[MAX_DIM];

        for(int i = 0; i < _dim; i++)
        {
            if(i != ax)
            {
                sub_shapes[i] = _shapes[i];
                sub_steps[i] = _steps[i];
            }
            else
            {
                sub_shapes[i] = end - start;
                sub_steps[i] = _steps[i] * step;
            }
        }
        return NDArray(_data + start * _steps[ax], _dim, sub_shapes, sub_steps, true, get_p_signature());
    }

    const size_t* get_p_signature() const
    {
        return _is_view? _p_signature: &_signature;
    }

    void _check_valid() const 
    {
        if(_is_view)
        {
            #ifdef DEBUG
            if(*get_p_signature() != get_signature() || !_data){
                panic("Array is invalid\n");
            }
            #endif
        }
    }

    NDArray<T> _contract() const
    {
        int new_shapes[MAX_DIM];
        int new_steps[MAX_DIM];
        int step = 0;
        int shape = 1;
        int new_dim = 0;

        for(int i = 0; i < _dim; i++)
        {
            step = _steps[i];
            shape *= _shapes[i];

            if(i == _dim-1 || _steps[i] != _steps[i+1] * _shapes[i+1])
            {
                new_shapes[new_dim] = shape;
                new_steps[new_dim]  = step;
                new_dim += 1;
                shape = 1;
                step = 0;
            }
        }
        
        return NDArray(_data, new_dim, new_shapes, new_steps, true, get_p_signature());
    }

    template<typename U>
    static void _contract_arrays(const NDArray& lhs, const NDArray<U>& rhs, NDArray& out_lhs, NDArray<U>& out_rhs)
    {
        int dim = std::max(lhs._dim, rhs._dim);
        int l_offset = dim - lhs._dim;
        int r_offset = dim - rhs._dim;

        int l_shapes[MAX_DIM];
        int r_shapes[MAX_DIM];
        int l_steps[MAX_DIM];
        int r_steps[MAX_DIM];

        for(int i = 0; i < lhs._dim; i++)
        {
            l_shapes[i+l_offset] = lhs._shapes[i];
            l_steps[i+l_offset] = lhs._steps[i];
        }

        for(int i = 0; i < rhs._dim; i++)
        {
            r_shapes[i+r_offset] = rhs._shapes[i];
            r_steps[i+r_offset] = rhs._steps[i];
        }

        int l_out_shapes[MAX_DIM];
        int r_out_shapes[MAX_DIM];
        int l_out_steps[MAX_DIM];
        int r_out_steps[MAX_DIM];
        int l_step = 0;
        int r_step = 0;
        int l_shape = 1;
        int r_shape = 1;
        int out_dim = 0;

        for(int i = 0; i < dim; i++)
        {
            l_step = l_steps[i];
            r_step = r_steps[i];
            l_shape *= l_shapes[i];
            r_shape *= r_shapes[i];

            bool split = false;
            if(i+1 == l_offset || i+1 == r_offset ||  i+1 == dim)
            {
                split = true;
            }

            if(i < dim-1 && i >= std::max(l_offset, r_offset) && 
                (l_shapes[i] != r_shapes[i] ||
                l_shapes[i+1] != r_shapes[i+1] || 
                l_steps[i] != l_steps[i+1] * l_shapes[i+1] ||
                r_steps[i] != r_steps[i+1] * r_shapes[i+1])
            )
            {
                split = true;
            }

            if(split)
            {
                l_out_shapes[out_dim] = l_shape;
                r_out_shapes[out_dim] = r_shape;
                l_out_steps[out_dim] = l_step;
                r_out_steps[out_dim] = r_step;
                out_dim += 1;
                l_shape = 1;
                r_shape = 1;
                l_step = 0;
                r_step = 0;
            }
        }
        
        l_offset = 0;
        while(l_offset < out_dim && l_out_shapes[l_offset] == 0) l_offset += 1;

        r_offset = 0;
        while(r_offset < out_dim && r_out_shapes[r_offset] == 0) r_offset += 1;
        
        out_lhs = NDArray(lhs._data, out_dim - l_offset,  l_out_shapes + l_offset, l_out_steps + l_offset, true, lhs.get_p_signature());
        out_rhs = NDArray<U>(rhs._data, out_dim - r_offset, r_out_shapes + r_offset, r_out_steps + r_offset, true, rhs.get_p_signature());
    }

    template<int func, typename U>
    void _unary_ops_compute(NDArray<U>& res) const
    {
        if(_dim == 1)
        {
            int size = _shapes[0];
            int step = _steps[0];
            int res_step = res._steps[0];
            for(int i = 0; i < size; i++)
            {
                res._data[i * res_step] = _unary_ops<func>(_data[i * step]);
            }
        }
        else
        {
            int min_ax = get_best_ax();
            for(int i = 0; i < _shapes[min_ax]; i++)
            {
                auto tmp = res._get_slice(min_ax, i);
                _get_slice(min_ax, i).template _unary_ops_compute<func>(tmp);
            }
        }
    }

    template<int ops, typename V, typename U>
    static void _bin_ops_compute(const NDArray& lhs, const V& rhs, NDArray<U>& res)
    {
        int dim = lhs.get_dim();

        if(dim == 1)
        {
            int size = lhs._shapes[0];
            int l_step = lhs._steps[0];
            int res_step = res._steps[0];

            for(int i = 0; i < size; i++)
            {
                res._data[i * res_step] = _bin_ops<ops>(lhs._data[i * l_step], rhs);
            }
            return;
        }

        int min_ax = lhs.get_best_ax();

        for(int i = 0; i < lhs._shapes[min_ax]; i++)
        {
            auto tmp = res._get_slice(min_ax, i);

            _bin_ops_compute<ops,V,U>(
                lhs._get_slice(min_ax, i),
                rhs,
                tmp
            );
        }
    }

    template<int ops, typename V, typename U>
    static void _bin_ops_compute(const V& rhs, const NDArray& lhs, NDArray<U>& res)
    {
        int dim = lhs.get_dim();

        if(dim == 1)
        {
            int size = lhs._shapes[0];
            int l_step = lhs._steps[0];
            int res_step = res._steps[0];

            for(int i = 0; i < size; i++)
            {
                res._data[i * res_step] = _bin_ops<ops>(rhs, lhs._data[i * l_step]);
            }
            return;
        }

        int min_ax = lhs.get_best_ax();

        for(int i = 0; i < lhs._shapes[min_ax]; i++)
        {
            auto tmp = res._get_slice(min_ax, i);

            _bin_ops_compute<ops,V,U>(
                rhs,
                lhs._get_slice(min_ax, i),
                tmp
            );
        }
    }

    template<int ops, typename V, typename U>
    static void _bin_ops_compute(const NDArray& lhs, const NDArray<V>& rhs, NDArray<U>& res)
    {
        int ldim = lhs.get_dim();
        int rdim = rhs.get_dim();
        int dim = std::max(ldim, rdim);

        if(dim == 1)
        {
            int l_size = lhs._shapes[0];
            int r_size = rhs._shapes[0];
            int l_step = lhs._steps[0];
            int r_step = rhs._steps[0];
            int res_step = res._steps[0];

            int size = std::max(l_size, r_size);

            if(l_size == 1)
            {
                for(int i = 0; i < size; i++)
                {
                    res._data[i * res_step] = _bin_ops<ops>(lhs._data[0], rhs._data[i * r_step]);
                }
                return;
            }
            else if(r_size == 1)
            {
                for(int i = 0; i < size; i++)
                {
                    res._data[i * res_step] = _bin_ops<ops>(lhs._data[i * l_step], rhs._data[0]);
                }
                return;
            }
            else
            {
                for(int i = 0; i < size; i++)
                {
                    res._data[i * res_step] = _bin_ops<ops>(lhs._data[i * l_step], rhs._data[i * r_step]);
                }
                return;
            }
        }

        if(dim > rdim)
        {
            int min_ax = lhs.get_min_ax(dim - rdim - 1);
            for(int i = 0; i < lhs._shapes[min_ax]; i++)
            {
                auto tmp = res._get_slice(min_ax, i);
                _bin_ops_compute<ops, V, U>(
                    lhs._get_slice(min_ax, i),
                    rhs,
                    tmp
                );
            }
            return;
        }

        if(dim > ldim)
        {
            int min_ax = rhs.get_min_ax(dim - ldim - 1);
            for(int i = 0; i < rhs._shapes[min_ax]; i++)
            {
                auto tmp = res._get_slice(min_ax, i);
                _bin_ops_compute<ops, V, U>(
                    lhs,
                    rhs._get_slice(min_ax, i),
                    tmp
                );
            }
            return;
        }

        int min_ax = get_best_ax(lhs, rhs);
        int l_shape = lhs._shapes[min_ax];
        int r_shape = rhs._shapes[min_ax];
        int shape = std::max(l_shape, r_shape);

        for(int i = 0; i < shape; i++)
        {
            int il = (i < l_shape) ? i : 0;
            int ir = (i < r_shape) ? i : 0;
            auto tmp = res._get_slice(min_ax, i);

            _bin_ops_compute<ops, V, U>(
                lhs._get_slice(min_ax, il),
                rhs._get_slice(min_ax, ir),
                tmp
            );
        }
    }

    template<typename U>
    bool _check_dim_assign(const NDArray<U>& oth) const
    {
        if(oth._dim > _dim) return false;
        int i1 = _dim - 1;
        int i2 = oth._dim - 1;
        while(i1 >= 0 && i2 >= 0)
        {
            if(_shapes[i1] != oth._shapes[i2] && oth._shapes[i2] > 1)
            {
                return false;
            }
            i1 -= 1;
            i2 -= 1;
        }
        return true;
    }

    template<typename U>
    bool _check_dim_bin_ops(const NDArray<U>& oth) const
    {
        int i1 = _dim - 1;
        int i2 = oth._dim - 1;
        while(i1 >= 0 && i2 >= 0)
        {
            if(_shapes[i1] != oth._shapes[i2] && _shapes[i1] > 1 && oth._shapes[i2] > 1)
            {
                return false;
            }
            i1 -= 1;
            i2 -= 1;
        }
        return true;
    }


    NDArray(const NDArray&) = delete;
    NDArray& operator=(const NDArray&) = delete;

    static int _calc_size(const int shapes[], int dim)
    {
        int size = 1;
        for(int i = 0; i < dim; i++) size *= shapes[i];
        return size;
    }

    static void _calc_steps(const int shapes[], int steps[], int dim)
    {
        int step = 1;
        for(int i = dim-1; i >= 0; i--)
        {
            steps[i] = step;
            step *= shapes[i];
        }
    }

    static NDArray _new(const int shapes[], int dim)
    {
        int size = _calc_size(shapes, dim);
        T* data = new T[size];
        int steps[MAX_DIM];
        _calc_steps(shapes, steps, dim);
        return NDArray(data, dim, shapes, steps);
    }

    static NDArray _full(const int shapes[], int dim, const T& val)
    {
        int size = _calc_size(shapes, dim);
        T* data = new T[size];
        for(int i = 0; i < size; i++) data[i] = val;
        int steps[MAX_DIM];
        _calc_steps(shapes, steps, dim);
        return NDArray(data, dim, shapes, steps);
    }

    NDArray _reshape(const int new_shapes[], int new_dim) const &
    {
        _check_valid();
        int min_ax = get_min_ax(_dim - 1);
        int new_steps[MAX_DIM];
        _calc_steps(new_shapes, new_steps, new_dim);
        int size = get_size();
        assert(new_steps[0] * new_shapes[0] == size);
        T* data = copy().into_raw();
        return NDArray(data, new_dim, new_shapes, new_steps);
    }

    NDArray _reshape2(const int new_shapes[], int new_dim)
    {
        _check_valid();
        int size = get_size();
        int new_steps[MAX_DIM];
        _calc_steps(new_shapes, new_steps, new_dim);
        assert(new_steps[0] * new_shapes[0] == size);

        if(!_is_view)
        {
            T* data = into_raw();
            return NDArray(data, new_dim, new_shapes, new_steps);
        }
        else
        {
            T* data = copy().into_raw();
            return NDArray(data, new_dim,  new_shapes, new_steps);
        }
    }

    void _calc_dim_expands(int new_shapes[], int ax)
    {
        new_shapes[ax] = 1;

        for(int i = 0; i < _dim; i++)
        {
            if(i < ax) new_shapes[i] = _shapes[i];
            else new_shapes[i+1] = _shapes[i];
        }
    }

    NDArray _view(const int new_shapes[], int new_dim) const
    {
        int new_size = _calc_size(new_shapes, new_dim);
        assert(!_is_view && new_size == get_size());
        int new_steps[MAX_DIM];
        _calc_steps(new_shapes, new_steps, new_dim);
        return NDArray(_data, new_dim, new_shapes, new_steps, true, get_p_signature());
    }

    static void _calc_bin_ops_shapes(const int l_shapes[], int l_dim, const int r_shapes[], int r_dim, int * out_shapes)
    {
        int dim = std::max(l_dim, r_dim);
        int l_offset = dim - l_dim;
        int r_offset = dim - r_dim;

        for(int i = 0; i < dim; i++)
        {
            int l_shape = i < l_offset ? 0 : l_shapes[i-l_offset];
            int r_shape = i < r_offset ? 0 : r_shapes[i-r_offset];
            out_shapes[i] = std::max(l_shape, r_shape);
        }
    }

    void _transpose_copy(NDArray<T>& result, const int* axes) const
    {
        if(_dim == 1)
        {
            result.copy_from(*this);
            return;
        }

        int min_ax = result.get_min_ax(_dim-1);
        int subaxes[MAX_DIM];
        for(int i = 0; i < _dim; i++)
        {
            int ax = axes[i];
            if(i < min_ax) subaxes[i] = axes[i];
            else if(i > min_ax) subaxes[i-1] = axes[i] - 1;
        }

        int ax0 = axes[min_ax];
        for(int i = 0; i < _shapes[ax0]; i++)
        {
            auto tmp = result._get_slice(min_ax, i);
            _get_slice(ax0, i)._transpose_copy(tmp, subaxes);
        }
    }
    
    NDArray _transpose(const int axes[]) const
    {
        _check_valid();

        for(int i = 0; i < _dim; i++)
        {
            boolean contained = false;
            for(int j = 0; j < _dim; j++)
            {
                if(axes[j] == i) 
                {
                    contained = true;
                    break;
                }
            }
            if(!contained)
            {
                panic("Axes missing index %d", i);
            }
        }

        int res_shapes[MAX_DIM];
        for(int i = 0; i < _dim; i++)
        {
            res_shapes[i] = _shapes[axes[i]];
        }
        auto result = NDArray::_new(res_shapes, _dim);
        _transpose_copy(result, axes);
        return result;
    }

    int get_scalar_index(const std::initializer_list<int>& indexes) const
    {
        assert(indexes.size() == get_dim());
        int idx = 0;
        int i = 0;
        for (int idx_i : indexes)
        {
            if(idx_i < 0) idx_i += _shapes[i];
            assert(idx_i >= 0 && idx_i < _shapes[i]);
            idx += idx_i * _steps[i];
            i += 1;
        }
        return idx;
    }

public:
    static NDArray empty()
    {
        return NDArray(NULL, 0, NULL, NULL);
    }
 
    NDArray()
    {

    }

    NDArray(const NDArray* oth): NDArray(oth->_data, oth->_shapes, oth->_steps, true, oth->get_p_signature())
    {
        
    }

    NDArray(NDArray&& oth)
    {
        _data = oth._data;
        _dim = oth._dim;
        for(int i = 0; i < _dim; i++)
        {
            _steps[i] = oth._steps[i];
            _shapes[i] = oth._shapes[i];
        }
        _is_view = oth._is_view;
        _p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
    }

    NDArray& operator=(NDArray&& oth)
    {
        if(_data && !_is_view)
        {
            delete[] _data;
        }
        _data = oth._data;
        _dim = oth._dim;
        for(int i = 0; i < _dim; i++)
        {
            _steps[i] = oth._steps[i];
            _shapes[i] = oth._shapes[i];
        }
        _is_view = oth._is_view;
        _p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
        return *this;
    }
    
    static NDArray full(const List<int>& shapes, const T& val)
    {
        return _full(shapes.get_data(), shapes.size(), val);
    }

    static NDArray ones(const List<int>& shapes)
    {
        return full(shapes, T(1));
    }

    static NDArray zeros(const List<int>& shapes)
    {
        return full(shapes, T(0));
    }

    NDArray copy() const
    {
        auto res = NDArray::_new(_shapes, _dim);
        res.copy_from(*this);
        return res;
    }

    NDArray view() const 
    {
        return NDArray(_data, _dim, _shapes, _steps, true, get_p_signature());
    }

    NDArray view(const List<int>& new_shapes) const
    {
        return _view(new_shapes.get_data(), new_shapes.size());
    }

    NDArray reshape(const List<int>& new_shapes) const &
    {
        _check_valid();
        return _reshape(new_shapes.get_data(), new_shapes.size());
    }

    NDArray reshape(const List<int>& new_shapes) &&
    {
        _check_valid();
        return _reshape2(new_shapes.get_data(), new_shapes.size());
    }

    NDArray flatten() const &
    {
        _check_valid();
        int size = get_size();
        int shapes[] = {size};
        return _reshape(shapes, 1);
    }

    NDArray flatten() &&
    {
        _check_valid();
        int size = get_size();
        int shapes[] = {size};
        return _reshape2(shapes, 1);
    }

    NDArray expand_dims(int ax) const &
    {
        _check_valid();
        if(ax < 0) ax += _dim;
        assert(ax >= 0 && ax < _dim);

        int new_shapes[MAX_DIM];
        _calc_dim_expands(new_shapes, ax);
        
        return _reshape(new_shapes, _dim+1);
    }

    NDArray expand_dims(int ax) &&
    {
        _check_valid();
        if(ax < 0) ax += _dim;
        assert(ax >= 0 && ax < _dim);

        int new_shapes[MAX_DIM];
        _calc_dim_expands(new_shapes, ax);
        
        return _reshape2(new_shapes, _dim+1);
    }

    T* get_data() const 
    {
        return _data;
    }

    T* into_raw()
    {
        _check_valid();
        assert(!_is_view);
        T* data = _data;
        _data = NULL;
        return data;
    }

    static NDArray from_raw(T* data, const List<int>& shapes)
    {
        int steps[MAX_DIM];
        _calc_steps(shapes.get_data(), steps, shapes.size());
        return NDArray(data, shapes.size(), shapes.get_data(), steps);
    }

    List<int> get_shapes() const
    {
        _check_valid();
        return List<int>(_shapes, _dim);
    }

    int get_shape(int idx) const 
    {
        _check_valid();
        if(idx < 0) idx += _dim;
        assert(idx >= 0 && idx < _dim);
        return _shapes[idx];
    }

    int get_dim() const 
    {
        _check_valid();
        return _dim;
    }

    int get_size() const
    {
        _check_valid();
        return _calc_size(_shapes, _dim);
    }

    template<int func>
    auto unary_ops() const
    {
        _check_valid();
        assert(func <= NUM_FUNC);
        using U = decltype(_unary_ops<func>(_data[0]));
        auto tmp = _contract();
        auto res = NDArray<U>::_new(tmp._shapes, tmp._dim);
        tmp.template _unary_ops_compute<func, U>(res);
        return std::move(res)._reshape(_shapes, _dim);
    }

    auto operator ~() const
    {
        return unary_ops<FUNC_NOT>();
    }

    auto operator -() const
    {
        return unary_ops<FUNC_NEG>();
    }

    auto abs() const 
    {
        return unary_ops<FUNC_ABS>();
    }

    auto round() const 
    {
        return unary_ops<FUNC_ROUND>();
    }

    auto ceil() const 
    {
        return unary_ops<FUNC_CEIL>();
    }

    auto floor() const 
    {
        return unary_ops<FUNC_FLOOR>();
    }

    auto square() const 
    {
        return unary_ops<FUNC_SQUARE>();
    }

    auto sqrt() const 
    {
        return unary_ops<FUNC_SQRT>();
    }

    auto exp() const 
    {
        return unary_ops<FUNC_EXP>();
    }

    auto log() const 
    {
        return unary_ops<FUNC_LOG>();
    }

    auto sin() const 
    {
        return unary_ops<FUNC_SIN>();
    }

    auto cos() const 
    {
        return unary_ops<FUNC_COS>();
    }

    auto tan() const 
    {
        return unary_ops<FUNC_TAN>();
    }

    auto asin() const 
    {
        return unary_ops<FUNC_ASIN>();
    }

    auto acos() const 
    {
        return unary_ops<FUNC_ACOS>();
    }

    auto atan() const 
    {
        return unary_ops<FUNC_ATAN>();
    }

    auto sinh() const 
    {
        return unary_ops<FUNC_SINH>();
    }

    auto cosh() const 
    {
        return unary_ops<FUNC_COSH>();
    }

    auto tanh() const 
    {
        return unary_ops<FUNC_TANH>();
    }

    auto asinh() const 
    {
        return unary_ops<FUNC_ASINH>();
    }

    auto acosh() const 
    {
        return unary_ops<FUNC_ACOSH>();
    }

    auto atanh() const 
    {
        return unary_ops<FUNC_ATANH>();
    }

    template<int ops, typename U, bool contract=true>
    void bin_ops_assign(const U& oth) requires std::copy_constructible<U>
    {
        _check_valid();

        auto tmp = contract ? _contract() : view();

        if(tmp.get_dim() == 1)
        {
            int size = tmp._shapes[0];
            int step = tmp._steps[0];
            for(int i = 0; i < size; i++)
            {
                _bin_ops_assign<ops>(tmp._data[i * step], oth);
            }
            return;
        }

        int min_ax = tmp.get_best_ax();

        for(int i = 0; i < tmp._shapes[min_ax]; i++)
        {
            tmp._get_slice(min_ax, i).template bin_ops_assign<ops, U, false>(oth);
        }
    }

    template<typename U>
    void set_all(const U& val) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_ASSIGN>(val);
    }

    template<typename U>
    void operator +=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_ADD>(oth);
    }

    template<typename U>
    void operator -=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_SUB>(oth);
    }

    template<typename U>
    void operator *=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_MUL>(oth);
    }

    template<typename U>
    void operator /=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_DIV>(oth);
    }

    template<typename U>
    void operator %=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_MOD>(oth);
    }

    template<typename U>
    void operator &=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_AND>(oth);
    }

    template<typename U>
    void operator |=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_OR>(oth);
    }

    template<typename U>
    void operator ^=(const U& oth) requires std::copy_constructible<U>
    {
        bin_ops_assign<OP_XOR>(oth);
    }

    template<int ops, typename U, bool contract=true>
    void bin_ops_assign(const NDArray<U>& oth)
    {
        assert(get_dim() >= oth.get_dim());
        if(contract && !_check_dim_assign(oth)) panic("Dimensions mismatch");

        auto lhs = view();
        auto rhs = oth.view();
        
        if(contract)
        {
            _contract_arrays(*this, oth, lhs, rhs);
        }

        if(lhs.get_dim() == 1)
        {
            int l_size = lhs.get_size();
            int r_size = rhs.get_size();
            if(r_size == 1)
            {
                for(int i = 0; i < l_size; i++)
                {
                    _bin_ops_assign<ops>(lhs._data[i * lhs._steps[0]], rhs._data[0]);
                }
            }
            else
            {
                int l_step = lhs._steps[0];
                int r_step = rhs._steps[0];
                for(int i = 0; i < l_size; i++)
                {
                    _bin_ops_assign<ops>(lhs._data[i * l_step], rhs._data[i * r_step]);
                }
            }
            return;
        }

        if(lhs.get_dim() > rhs.get_dim())
        {
            int min_ax = lhs.get_min_ax(lhs.get_dim() - rhs.get_dim() - 1);
            for(int i = 0; i < lhs._shapes[min_ax]; i++)
            {
                lhs._get_slice(min_ax, i).template bin_ops_assign<ops, U, false>(rhs);
            }
            return;
        }

        int min_ax =  get_best_ax(lhs, rhs);

        for(int i = 0; i < lhs._shapes[min_ax]; i++)
        {
            int ir = (i < rhs._shapes[min_ax]) ? i : 0;
            lhs._get_slice(min_ax, i).template bin_ops_assign<ops, U, false>(rhs._get_slice(min_ax, ir));
        }
    }


    template<typename U>
    void copy_from(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_ASSIGN, U>(oth);
    }

    template<typename U>
    NDArray<U> as_type()
    {
        auto res = NDArray<U>::_new(_shapes, _dim);
        res.copy_from(*this);
        return res;
    }

    template<typename U>
    void operator +=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_ADD, U>(oth);
    }

    template<typename U>
    void operator -=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_SUB, U>(oth);
    }

    template<typename U>
    void operator *=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_MUL, U>(oth);
    }

    template<typename U>
    void operator /=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_DIV, U>(oth);
    }

    template<typename U>
    void operator %=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_MOD, U>(oth);
    }

    template<typename U>
    void operator &=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_AND, U>(oth);
    }

    template<typename U>
    void operator |=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_OR, U>(oth);
    }

    template<typename U>
    void operator ^=(const NDArray<U>& oth)
    {
        bin_ops_assign<OP_XOR, U>(oth);
    }

    template<int ops, typename V>
    auto bin_ops(const V& rhs) const requires std::copy_constructible<V>
    {
        _check_valid();
        assert(ops <= NUM_OP);
        auto tmp = _contract();
        using U = decltype(_bin_ops<ops>(_data[0], rhs));
        auto res = NDArray<U>::_new(tmp._shapes, tmp._dim);
        _bin_ops_compute<ops, V, U>(tmp, rhs, res);
        return std::move(res)._reshape(_shapes, _dim);
    }

    template<int ops, typename V>
    auto bin_ops_reverse(const V& rhs) const requires std::copy_constructible<V>
    {
        _check_valid();
        assert(ops <= NUM_OP);
        auto tmp = _contract();
        using U = decltype(_bin_ops<ops>(_data[0], rhs));
        auto res = NDArray<U>::_new(tmp._shapes, tmp._dim);
        _bin_ops_compute<ops, V, U>(rhs, tmp, res);
        return std::move(res)._reshape(_shapes, _dim);
    }

    template<typename  V>
    auto minimum(const V& rhs) const requires std::copy_constructible<V>
    {
        return bin_ops<OP_MINIMUM, V>(rhs);
    }

    template<typename  V>
    auto maximum(const V& rhs) const requires std::copy_constructible<V>
    {
        return bin_ops<OP_MAXIMUM, V>(rhs);
    }

    template<int ops, typename V>
    auto bin_ops(const NDArray<V>& oth) const
    {
        _check_valid();

        if(!_check_dim_bin_ops(oth)) panic("Dimensions mismatch");

        assert(ops <= NUM_OP);

        using U = decltype(_bin_ops<ops>(_data[0], oth._data[0]));
        
        NDArray<T> lhs;
        NDArray<V> rhs;
        _contract_arrays(*this, oth, lhs, rhs);

        int res_shapes[MAX_DIM];
        _calc_bin_ops_shapes(lhs._shapes, lhs._dim, rhs._shapes, rhs._dim, res_shapes);
        auto res = NDArray<U>::_new(res_shapes, std::max(lhs._dim, rhs._dim));
        _bin_ops_compute<ops, V, U>(lhs, rhs, res);

        _calc_bin_ops_shapes(_shapes, _dim, oth._shapes, oth._dim, res_shapes);
        return std::move(res)._reshape2(res_shapes, std::max(_dim, oth._dim));
    }

    template<typename  V>
    auto operator +(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_ADD, V>(rhs);
    }

    template<typename  V>
    auto operator -(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_SUB, V>(rhs);
    }

    template<typename  V>
    auto operator *(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_MUL, V>(rhs);
    }

    template<typename  V>
    auto operator /(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_DIV, V>(rhs);
    }

    template<typename  V>
    auto operator %(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_MOD, V>(rhs);
    }

    template<typename  V>
    auto pow(const NDArray<V>& rhs) const 
    {
        return bin_ops<OP_POW, V>(rhs); 
    }

    template<typename  V>
    auto operator &(const NDArray<V>& rhs)  const
    {
        return bin_ops<OP_AND, V>(rhs);
    }

    template<typename  V>
    auto operator |(const NDArray<V>& rhs) const
    {
        return bin_ops<T, OP_OR>(rhs);
    }

    template<typename  V>
    auto operator ^(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_XOR, V>(rhs);
    }

    template<typename  V>
    auto operator >(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_GT, V>(rhs);
    }

    template<typename  V>
    auto operator >=(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_GTE, V>(rhs);
    }

    template<typename  V>
    auto operator <(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_LT, V>(rhs);
    }

    template<typename  V>
    auto operator <=(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_LTE, V>(rhs);
    }

    template<typename  V>
    auto operator ==(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_EQ, V>(rhs);
    }

    template<typename  V>
    auto operator !=(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_NE, V>(rhs);
    }

    template<typename  V>
    auto minimum(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_MINIMUM, V>(rhs);
    }

    template<typename  V>
    auto maximum(const NDArray<V>& rhs) const
    {
        return bin_ops<OP_MAXIMUM, V>(rhs);
    }

    double percentile(double p) const
    {
        int size = get_size();
        assert(size > 0);
        auto arr = copy();
        arr.sort();
        double index = (p/100) * (size-1);
        int index_int = (int) index;
        double delta = index - index_int;
        if(index_int == size-1)
        {
            return (double) _data[index_int];
        }
        return (1.0 - delta) * arr._data[index_int] + delta * arr._data[index_int+1];
    }

    template<int ops>
    auto reduce_ops() const -> decltype(_reduce_ops<ops>(_data[0], _data[0]))
    {
        if constexpr(ops == REDUCE_ARGMIN)
        {
            return argmin();
        }
        else if constexpr(ops == REDUCE_ARGMAX)
        {
            return argmax();
        }
        else if constexpr(ops == REDUCE_MEAN)
        {
            return reduce_mean();
        }
        else if constexpr(ops == REDUCE_MEDIAN)
        {
            return reduce_median();
        }
        else if constexpr(ops == REDUCE_STD)
        {
            return reduce_std();
        }
        else
        {
            using U = decltype(_reduce_ops<ops>(_data[0], _data[0]));
            auto acc = U(0);
            if constexpr(ops == REDUCE_MIN || ops == REDUCE_MAX)
            {
                assert(_shapes[0] > 0);
                acc = _data[0];
            }
            if constexpr(ops == REDUCE_AND) acc = U(1);

            auto tmp = _contract();

            if(tmp.get_dim() == 1)
            {
                for(int i = 0; i < tmp._shapes[0]; i++)
                {
                    acc = _reduce_ops<ops>(acc, tmp._data[i * tmp._steps[0]]);
                }
            }
            else
            {
                int min_ax = tmp.get_best_ax();
                for(int i = 0; i < tmp._shapes[min_ax]; i++)
                {
                    U tmp_acc = tmp._get_slice(min_ax, i).template reduce_ops<ops>();
                    acc = _reduce_ops<ops>(
                        acc,
                        tmp_acc
                    );
                }
            }
            return acc;
        }
    }

    auto reduce_max() const 
    {
        return reduce_ops<REDUCE_MAX>();
    }

    T reduce_min() const 
    {
        return reduce_ops<REDUCE_MIN>();
    }

    auto reduce_sum() const 
    {
        return reduce_ops<REDUCE_SUM>();
    }

    double reduce_mean() const 
    {
        int size = get_size();
        assert(size > 0);
        return (double) reduce_sum() / size;
    }

    double reduce_median() const
    {
        assert(get_size() > 0);
        return percentile(50);
    }

    double reduce_std() const 
    {
        int size = get_size();
        assert(size > 0);

        auto arr_square = square();
        double m2 = arr_square.reduce_mean();
        double m = reduce_mean();
        return ::sqrt(m2 - m*m);
    }

    auto reduce_and() const 
    {
        return reduce_ops<REDUCE_AND>();
    }

    auto reduce_or() const 
    {
        return reduce_ops<REDUCE_OR>();
    }

    auto reduce_xor() const 
    {
        return reduce_ops<REDUCE_XOR>();
    }

    int argmax() const 
    {
        if(get_dim() > 1)
        {
            return flatten().argmax();
        }
        int size = _shapes[0];

        if(size == 0)
        {
            panic("Max of empty array\n");
        }

        T max_val = _data[0];
        int imax = 0;

        for(int i = 1; i < size; i++)
        {
            T tmp = _data[i * _steps[0]];
            if(tmp > max_val)
            {
                imax = i;
                max_val = tmp;
            }
        }
        return imax;
    }

    int argmin() const 
    {
        if(get_dim() > 1)
        {
            return flatten().argmin();
        }
        int size = _shapes[0];

        if(size == 0)
        {
            panic("Max of empty array\n");
        }

        T min_val = _data[0];
        int imin = 0;

        for(int i = 1; i < size; i++)
        {
            T tmp = _data[i * _steps[0]];
            if(tmp < min_val)
            {
                imin = i;
                min_val = tmp;
            }
        }
        return imin;
    }

    template<int ops, typename U>
    void _get_reduce_result(NDArray<U>& result, int axis) const
    {   
        int dim = get_dim();
        assert(dim > 1);
        int best_ax = -1;
        for(int i = 0; i < dim; i++)
        {
            if(i == axis) continue;
            if(best_ax < 0 || _shapes[i] < _shapes[best_ax])
            {
                best_ax = i;
            }
        }

        if(dim == 2)
        {
            for(int i = 0; i < _shapes[best_ax]; i++)
            {
                result.at(i) = _get_slice(best_ax, i).template reduce_ops<ops>();
            }
        }
        else 
        {
            int n = _shapes[best_ax];
            int res_ax = best_ax;
            if(best_ax < axis) axis -= 1;else res_ax -= 1;
            for(int i = 0; i < n; i++)
            {
                auto tmp = result._get_slice(res_ax, i);
                _get_slice(best_ax, i).template _get_reduce_result<ops, U>(tmp, axis);
            }
        }
    }

    template<int ops>
    auto reduce_ops(int axis) const -> NDArray<decltype(_reduce_ops<ops>(_data[0], _data[0]))>
    {
        _check_valid();
        if(axis < 0) axis += _dim;
        assert(_dim > 1 && axis >= 0 && axis < _dim);

        int subshapes[MAX_DIM];

        for(int i = 0; i < _dim; i++)
        {
            if(i < axis) subshapes[i] = _shapes[i];
            else if(i > axis) subshapes[i-1] = _shapes[i];
        }

        using U = decltype(_reduce_ops<ops>(_data[0], _data[0]));
        auto result = NDArray<U>::_new(subshapes, _dim-1);
        _get_reduce_result<ops, U>(result, axis);
        return result;
    }


    auto argmax(int axis) const
    {
        return reduce_ops<REDUCE_ARGMAX>(axis);
    }

    auto argmin(int axis) const
    {
        return reduce_ops<REDUCE_ARGMIN>(axis);
    }

    auto reduce_sum(int axis) const
    {
        return reduce_ops<REDUCE_SUM>(axis);
    } 

    auto reduce_min(int axis) const
    {
        return reduce_ops<REDUCE_MIN>(axis);
    }

    auto reduce_max(int axis) const
    {
        return reduce_ops<REDUCE_MAX>(axis);
    }

    auto reduce_mean(int axis) const
    {
        return reduce_ops<REDUCE_MEAN>(axis);
    } 

    auto reduce_median(int axis) const
    {
        return reduce_ops<REDUCE_MEDIAN>(axis);
    }

    auto reduce_std(int axis) const
    {
        return reduce_ops<REDUCE_STD>(axis);
    }

    auto reduce_and(int axis) const
    {
        return reduce_ops<REDUCE_AND>(axis);
    }

    auto reduce_or(int axis) const
    {
        return reduce_ops<REDUCE_OR>(axis);
    } 


    void sort(int axis=-1)
    {
        _check_valid();
        int dim = get_dim();
        if(axis < 0) axis += dim;
        assert(axis >= 0 && axis < dim);

        if(dim > 1)
        {
            int best_ax = -1;
            for(int i = 0; i < dim; i++)
            {
                if(i == axis) continue;
                if(best_ax < 0 || _shapes[i] < _shapes[best_ax])
                {
                    best_ax = i;
                }
            }
            if(best_ax < axis) axis -= 1;
            for(int i = 0; i < _shapes[best_ax]; i++)
            {
                _get_slice(best_ax, i).sort(axis);
            }
        }
        else if(_shapes[0] > 0)
        {
            quick_sort(_data, _steps[0], NULL, 0, _shapes[0] - 1);
        }
    }

    NDArray<int> argsort(int axis=-1) const
    {
        int dim = get_dim();
        if(axis < 0) axis += dim;
        assert(axis >= 0 && axis < dim);
        auto result = NDArray<int>::_new(_shapes, _dim);

        if(dim > 1)
        {
            int best_ax = -1;
            for(int i = 0; i < dim; i++)
            {
                if(i == axis) continue;
                if(best_ax < 0 || _shapes[i] < _shapes[best_ax])
                {
                    best_ax = i;
                }
            }
            if(best_ax < axis) axis -= 1;
            for(int i = 0; i < _shapes[best_ax]; i++)
            {
                auto tmp = result._get_slice(best_ax, i);
                tmp.copy_from(_get_slice(best_ax, i).argsort(axis));
            }
        }
        else
        {
            for(int i=0; i < _shapes[0]; i++) result._data[i] = i;
            auto tmp = copy();
            quick_sort(tmp._data, 1, result._data, 0, _shapes[0]-1);
            
        }
        return result;
    }

    T& at(int idx) 
    {
        assert(get_dim() == 1);
        if(idx < 0) idx += _shapes[0];
        assert(idx >= 0 && idx < _shapes[0]);
        return _data[idx * _steps[0]];
    }

    const T& at(int idx) const
    {
        assert(get_dim() == 1);
        if(idx < 0) idx += _shapes[0];
        assert(idx >= 0 && idx < _shapes[0]);
        return _data[idx * _steps[0]];
    }

    T& at(const std::initializer_list<int>& indexes)
    {
        return _data[get_scalar_index(indexes)];
    }

    const T& at(const std::initializer_list<int>& indexes) const
    {
        return _data[get_scalar_index(indexes)];
    }

    /*
    template<typename ...Args>
    NDArray operator[](const Args&... args)
    {
        return at({args...});
    }
    */
    T& operator[](const int& idx) 
    {
        return at(idx);
    }

    const T& operator[](const int& idx) const
    {
        return at(idx);
    }

    NDArray operator[](const Index& index) const
    {
        return slice(&index, 1);
    }
    
    NDArray operator[](const List<Index>& indexes) const
    {
        return slice(indexes.get_data(), indexes.size());
    }

    NDArray slice(const Index* indexes, int n) const
    {
        int n_list = 0;
        for(int i = 0; i < n; i++)
        {
            if(indexes[i].is_list())
            {
                n_list += 1;
            }
        }
        if(n_list >= 2) return slice3(indexes, n);
        if(n_list == 1) return slice2(indexes, n);
        return slice1(indexes, n);
    }

    NDArray slice1(const Index* indexes, int n) const
    {
        int sub_shapes[MAX_DIM];
        int sub_steps[MAX_DIM];
        int sub_dim = 0;

        int offset = 0;
        int ax = 0;

        for(int i = 0; i < n; i++)
        {
            if(indexes[i].is_none())
            {
                sub_shapes[sub_dim] = 1;
                sub_steps[sub_dim] = 0;
                sub_dim += 1;
            }
            else
            {
                assert(ax < _dim);

                if(indexes[i].is_scalar())
                {
                    int idx_i = indexes[i].get_scalar();
                    if(idx_i < 0) idx_i += _shapes[ax];
                    assert(idx_i >= 0 && idx_i < _shapes[ax]);
                    offset += idx_i * _steps[ax];    
                }
                else if(indexes[i].is_slice())
                {
                    Slice index = *indexes[i].get_slice();
                    index.resolve_index(_shapes[ax]);
                    offset += index.start.value * _steps[ax];
                    sub_shapes[sub_dim] = (index.stop.value - index.start.value) / index.step.value;

                    if(sub_dim > 0 && sub_steps[sub_dim-1] == 0) 
                        sub_steps[sub_dim-1] = _steps[ax] * index.step.value;

                    sub_steps[sub_dim] = _steps[ax] * index.step.value;
                    sub_dim += 1;
                }
                ax += 1;
            }
        }

        for(; ax < _dim; ax++)
        {
            sub_shapes[sub_dim] = _shapes[ax];
            if(sub_dim > 0 && sub_steps[sub_dim-1] == 0) 
                sub_steps[sub_dim-1] = _steps[ax];

            sub_steps[sub_dim] = _steps[ax];
            sub_dim += 1;
        }

        return NDArray(_data + offset, sub_dim, sub_shapes, sub_steps, true, get_p_signature());
    }

    void _copy_slice(const Index* indexes, int index_size, NDArray& res) const
    {
        int min_ax = get_best_ax(indexes, index_size);

        Index sub_indexes[MAX_DIM];
        for(int i = 0; i < index_size; i++)
        {
            if(i < min_ax) sub_indexes[i] = indexes[i];
            else if(i > min_ax) sub_indexes[i-1] = indexes[i];
        }

        if(indexes[min_ax].is_list())
        {
            assert(indexes[min_ax].get_max_idx() < 0 || indexes[min_ax].get_max_idx() == _shapes[min_ax]);
            auto idx = indexes[min_ax].get_indexes();
            int idx_len = indexes[min_ax].get_index_length();

            if(index_size > 1)
            {
                for(int i = 0; i < idx_len; i++)
                {
                    auto tmp = res._get_slice(min_ax, i);
                    _get_slice(min_ax, idx[i])._copy_slice(sub_indexes, index_size-1, tmp);
                }
            }
            else
            {
                for(int i = 0; i < idx_len; i++)
                {
                    res._data[i] = _data[idx[i] * _steps[0]];
                }
            }
        }
        else
        {
            auto slice = indexes[min_ax].get_slice();
            int start = slice->start.value;
            int stop = slice->stop.value;
            int step = slice->step.value;
            if(index_size > 1)
            {
                int i = 0;
                while(start != stop)
                {
                    auto tmp = res._get_slice(min_ax, i);
                    _get_slice(min_ax, start)._copy_slice(sub_indexes, index_size-1, tmp);
                    i += 1;
                    start += step;
                }
            }
            else
            {
                int i = 0;
                while(start != stop)
                {
                    res._data[i] = _data[start * _steps[0]];
                    i += 1;
                    start += step;
                }
            }
        }
    }

    NDArray slice2(const Index* indexes, int n) const
    {
        int output_shapes[MAX_DIM];
        int tmp_shapes[MAX_DIM];
        Index indexes1[MAX_DIM];
        Index indexes2[MAX_DIM];
        int output_dim = 0, tmp_dim = 0, indexes1_size = 0, indexes2_size = 0;
        int ax = 0;

        for(int i = 0; i < n; i++)
        {
            Index index = indexes[i];
            if(index.is_none())
            {
                output_shapes[output_dim] = 1;
                output_dim += 1;
            }
            else
            {
                assert(ax < _dim);
                if(index.is_slice())
                {
                    Slice* slice = index.get_slice();
                    slice->resolve_index(_shapes[ax]);
                    int sz = (slice->stop.value - slice->start.value) / (slice->step.value);
                    output_shapes[output_dim] = sz;
                    tmp_shapes[tmp_dim] = sz;
                    output_dim += 1;
                    tmp_dim += 1;
                }
                else if(index.is_list())
                {
                    int sz = index.get_index_length();
                    output_shapes[output_dim] = sz;
                    tmp_shapes[tmp_dim] = sz;
                    output_dim += 1;
                    tmp_dim += 1;
                }
             
                if(index.is_scalar())
                {
                    indexes1[indexes1_size] = index;
                    indexes2[indexes2_size] = Slice::all();
                    indexes1_size += 1;
                    indexes2_size += 1;
                }
                else
                {
                    indexes1[indexes1_size] = Slice::all();
                    indexes2[indexes2_size] = index;
                    indexes1_size += 1;
                    indexes2_size += 1;
                }
                ax += 1;
            }
        }

        for(; ax < _dim; ax++)
        {
            output_shapes[output_dim] = _shapes[ax];
            tmp_shapes[tmp_dim] = _shapes[ax];
            output_dim += 1;
            tmp_dim += 1;
        }

        auto res = NDArray::_new(output_shapes, output_dim);
        auto res_tmp = res._view(tmp_shapes, tmp_dim);

        auto tmp = slice(indexes1, indexes1_size);
        tmp._copy_slice(indexes2, indexes2_size, res_tmp);
        return res;
    }

    NDArray slice3(const Index* indexes, int n) const
    {
        int output_shapes[MAX_DIM];
        int output_dim = 1;
        output_shapes[0] = 0;
        int ax = 0;

        for(int i = 0; i < n; i++)
        {
            Index index = indexes[i];
            if(index.is_none())
            {
                output_shapes[output_dim] = 1;
                output_dim += 1;
            }
            else
            {
                assert(ax < _dim);
                if(index.is_slice())
                {
                    Slice* slice = index.get_slice();
                    slice->resolve_index(_shapes[ax]);
                    int sz = (slice->stop.value - slice->start.value) / (slice->step.value);
                    output_shapes[output_dim] = sz;
                    output_dim += 1;
                }
             
                ax += 1;
            }
        }

        for(; ax < _dim; ax++)
        {
            output_shapes[output_dim] = _shapes[ax];
            output_dim += 1;
        }

        int n_elem = -1;
        int ifirst = -1;
        for(int i = 0; i < n; i++)
        {
            if(indexes[i].is_list()) 
            {
                int sz = indexes[i].get_index_length();
                if(sz == 0)
                {
                    n_elem = 0;
                    break;
                }

                if(n_elem < 0)
                {
                    n_elem = sz;
                    ifirst = i;
                }
                
                if(sz != 1 && sz != n_elem)
                {
                    panic("Indexes in axis 0 and %d do not have same number of elements: %d != %d", i, n_elem, sz);
                }

                int max_idx = indexes[i].get_max_idx();
                if(max_idx >= 0 && max_idx != _shapes[i])
                {
                    panic("Mask size and array size in axis %d not equal: %d != %d", i, indexes[i].get_max_idx(), _shapes[i]);
                }
            }
        }

        int dim = get_dim();

        output_shapes[0] = n_elem;
        auto res = NDArray::_new(output_shapes, output_dim);
        Index subindexes[MAX_DIM];
        for(int i = 0; i < n; i++) subindexes[i] = indexes[i];

        for(int k = 0; k < n_elem; k++)
        {
            for(int i = 0; i < n; i++)
            {
                if(indexes[i].is_list()) 
                {
                    int* idx_i = indexes[i].get_indexes();
                    int sz = indexes[i].get_index_length();
                    subindexes[i] = Index(idx_i[k < sz ? k : 0]);
                }
            }
            res._get_slice(0, k).copy_from(slice(subindexes, n));
        }
        return res;
    }

    NDArray transpose(const List<int>& axes) const
    {
        assert(axes.size() == _dim);
        return transpose(axes.get_data(), axes.size());
    }

    static NDArray concatenate(const List<NDArray>& lst, int ax=0)
    {
        if(lst.size() == 0) return empty();
        int dim = lst[0].get_dim();
        if(ax < 0) ax += dim;
        assert(ax >= 0 && ax < dim);

        for(int i = 1; i < lst.size(); i++)
        {
            if(lst[i].get_dim() != dim)
            {
                panic("Cannot concatenate: Array 0 has %d dimension, but array %d has %d dimension", i, dim, lst[i].get_dim());
            }
        }

        for(int d = 0; d < dim; d++)
        {
            int sh = lst[0]._shapes[d];
            if(d == ax) continue;
            for(int i = 1; i < lst.size(); i++)
            {
                if(lst[i]._shapes[d] != sh)
                {
                    panic("Cannot concatenate: Array sizes accross dimension %d mismatch: array 0 has %d elements, but array %d has %d elements",
                         d, sh, i, lst[i]._shapes[d]);
                }
            }
        }
        
        int concat_size = 0;
        for(auto& arr : lst) concat_size += arr._shapes[ax];
        int concat_shapes[MAX_DIM];

        for(int d = 0; d < dim; d++)
        {
            if(d == ax) 
            {
                concat_shapes[d] = concat_size;
            }
            else
            {
                concat_shapes[d] = lst[0]._shapes[d];
            }
        }
        auto result = NDArray::_new(concat_shapes, dim);

        int offset = 0;
        for(auto& arr: lst)
        {
            auto tmp = result._get_slice(ax, offset, offset + arr._shapes[ax]);
            tmp.copy_from(arr);
            offset += arr._shapes[ax];
        }

        return result;
    }

    static NDArray stack(const List<NDArray>& lst, int ax=0)
    {
        if(lst.size() == 0) return empty();
        int dim = lst[0].get_dim();
        if(ax < 0) ax += dim+1;
        assert(ax >= 0 && ax <= dim);

        for(int i = 1; i < lst.size(); i++)
        {
            if(lst[i].get_dim() != dim)
            {
                panic("Cannot stack: Array 0 has %d dimension, but array %d has %d dimension", i, dim, lst[i].get_dim());
            }
        }

        for(int d = 0; d < dim; d++)
        {
            int sh = lst[0]._shapes[d];
            for(int i = 1; i < lst.size(); i++)
            {
                if(lst[i]._shapes[d] != sh)
                {
                    panic("Cannot stack: Array sizes accross dimension %d mismatch: array 0 has %d elements, but array %d has %d elements",
                         d, sh, i, lst[i]._shapes[d]);
                }
            }
        }
        
        int stack_shapes[MAX_DIM];
        stack_shapes[ax] = lst.size();
        for(int d = 0; d < dim; d++)
        {
            if(d < ax) stack_shapes[d] = lst[0]._shapes[d];
            else stack_shapes[d+1] = lst[0]._shapes[d];
        }

        auto result = NDArray::_new(stack_shapes, dim+1);

        for(int i = 0; i < lst.size(); i++)
        {
            auto tmp = result._get_slice(ax, i);
            tmp.copy_from(lst[i]);
        }

        return result;
    }

    ~NDArray()
    {
        _signature = 0;
        if(!_is_view && _data) {
            delete[] _data;
        }
    }
};

template<typename T>
static std::ostream & operator<<(std::ostream & os, const NDArray<T>& arr) 
{
    int sz = arr.get_shape(0);
    os << "[";
    if(arr.get_dim() == 1)
    {
        for(int i = 0; i < sz; i++)
        {
            os << arr.at(i);
            if(i + 1 < sz) os << ", ";
        }    
    }
    else
    {
        for(int i = 0; i < sz; i++)
        {
            os << arr.get_slice(0, i);
            if(i+1 < sz) os << ", ";
        }
    }
    os << "]";
    
    return os;
}

template<typename T, typename  V>
auto operator +(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_ADD, V>(rhs);
}

template<typename T, typename  V>
auto operator -(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_SUB, V>(rhs);
}

template<typename T, typename  V>
auto operator *(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_MUL, V>(rhs);
}

template<typename T, typename  V>
auto operator /(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_DIV, V>(rhs);
}

template<typename T, typename  V>
auto operator %(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_MOD, V>(rhs);
}

template<typename T, typename  V>
auto pow(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_POW, V>(rhs); 
}

template<typename T, typename  V>
auto operator &(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_AND, V>(rhs);
}

template<typename T, typename  V>
auto operator |(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<T, OP_OR>(rhs);
}

template<typename T, typename  V>
auto operator ^(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_XOR, V>(rhs);
}

template<typename T, typename  V>
auto operator >(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_GT, V>(rhs);
}

template<typename T, typename  V>
auto operator >=(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_GTE, V>(rhs);
}

template<typename T, typename  V>
auto operator <(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_LT, V>(rhs);
}

template<typename T, typename  V>
auto operator <=(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_LTE, V>(rhs);
}

template<typename T, typename  V>
auto operator ==(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_EQ, V>(rhs);
}

template<typename T, typename  V>
auto operator !=(const NDArray<T>& arr, const V& rhs) requires std::copy_constructible<V>
{
    return arr.template bin_ops<OP_NE, V>(rhs);
}

//  Binops reverse
template<typename T, typename  V>
auto operator +(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_ADD, V>(rhs);
}

template<typename T, typename  V>
auto operator -(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_SUB, V>(rhs);
}

template<typename T, typename  V>
auto operator *(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_MUL, V>(rhs);
}

template<typename T, typename  V>
auto operator /(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_DIV, V>(rhs);
}

template<typename T, typename  V>
auto operator %(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_MOD, V>(rhs);
}

template<typename T, typename  V>
auto pow(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_POW, V>(rhs); 
}

template<typename T, typename  V>
auto operator &(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_AND, V>(rhs);
}

template<typename T, typename  V>
auto operator |(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<T, OP_OR>(rhs);
}

template<typename T, typename  V>
auto operator ^(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_XOR, V>(rhs);
}

template<typename T, typename  V>
auto operator >(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_GT, V>(rhs);
}

template<typename T, typename  V>
auto operator >=(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_GTE, V>(rhs);
}

template<typename T, typename  V>
auto operator <(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_LT, V>(rhs);
}

template<typename T, typename  V>
auto operator <=(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_LTE, V>(rhs);
}

template<typename T, typename  V>
auto operator ==(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_EQ, V>(rhs);
}

template<typename T, typename  V>
auto operator !=(const V& rhs, const NDArray<T>& arr) requires std::copy_constructible<V>
{
    return arr.template bin_ops_reverse<OP_NE, V>(rhs);
}

template<typename T>
class Array : public NDArray<T>
{
public:
    Array() : NDArray<T>()
    {
        
    }

    Array(NDArray<T>&& oth)
    {
        assert(oth._dim == 1);
        this->_data = oth._data;
        this->_dim = 1;
        this->_steps[0] = oth._steps[0];
        this->_shapes[0] = oth._shapes[0];
        this->_is_view = oth._is_view;
        this->_p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
    }

    Array(const List<T>& lst)
    {
        int size = lst.size();
        T* ptr = new T[size];
        const T* data = lst.get_data();
        for(int i = 0; i < size; i++) ptr[i] = data[i];
        this->_data = ptr;
        this->_dim = 1;
        this->_steps[0] = 1;
        this->_shapes[0] = size;
        this->_is_view = false;
    }
};

template<typename T>
class Array2 : public NDArray<T>
{
public:
    Array2() : NDArray<T>()
    {

    }
    Array2(NDArray<T>&& oth)
    {
        assert(oth.get_dim() == 2);
        this->_data = oth._data;
        this->_dim = 2;
        this->_steps[0] = oth._steps[0];
        this->_steps[1] = oth._steps[1];
        this->_shapes[0] = oth._shapes[0];
        this->_shapes[1] = oth._shapes[1];
        this->_is_view = oth._is_view;
        this->_p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
    }
};

template<typename T>
class Array3 : public NDArray<T>
{
public:
    Array3() : NDArray<T>()
    {
        
    }

    T& at(int i, int j, int k)
    {
        return this->_data[i * this->_steps[0] + j * this->_steps[1] + k * this->_steps[2]];
    }

    const T& at(int i, int j, int k) const
    {
        return this->_data[i * this->_steps[0] + j * this->_steps[1] + k * this->_steps[2]];
    }

    Array3(NDArray<T>&& oth)
    {
        assert(oth.get_dim() == 3);
        this->_data = oth._data;
        this->_dim = 3;
        this->_steps[0] = oth._steps[0];
        this->_steps[1] = oth._steps[1];
        this->_steps[2] = oth._steps[2];
        this->_shapes[0] = oth._shapes[0];
        this->_shapes[1] = oth._shapes[1];
        this->_shapes[2] = oth._shapes[2];
        this->_is_view = oth._is_view;
        this->_p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
    }
};

template<typename T>
class Array4 : public NDArray<T>
{
public:
    Array4() : NDArray<T>()
    {
        
    }

    Array4(NDArray<T>&& oth)
    {
        assert(oth.get_dim() == 4);
        this->_data = oth._data;
        this->_dim = 4;
        this->_steps[0] = oth._steps[0];
        this->_steps[1] = oth._steps[1];
        this->_steps[2] = oth._steps[2];
        this->_steps[3] = oth._steps[3];
        this->_shapes[0] = oth._shapes[0];
        this->_shapes[1] = oth._shapes[1];
        this->_shapes[2] = oth._shapes[2];
        this->_shapes[3] = oth._shapes[3];
        this->_is_view = oth._is_view;
        this->_p_signature = oth._p_signature;
        this->_signature = oth._signature;
        oth._data = NULL;
    }
};
#endif