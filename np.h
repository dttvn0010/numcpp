#ifndef __NP_H__
#define __NP_H__
#include "base.h"
#include "nd_array.h"
#include <cassert>

namespace np 
{
    template<typename T>
    NDArray<T> full(const List<int>& shapes, const T& val)
    {
        return NDArray<T>::full(shapes, val);
    }

    template<typename T>
    NDArray<T> zeros(const List<int>& shapes)
    {
        return full(shapes, T(0));
    }

    template<typename T>
    NDArray<T> ones(const List<int>& shapes)
    {
        return full(shapes, T(1));
    }

    template<typename  T>
    NDArray<T> array(const List<T>& lst, const List<int>& shapes)
    {
        int size = 1;
        for(int sh : shapes) size *= sh;
        assert(size == lst.size());
        T* ptr = new T[size];
        const T* data = lst.get_data();
        for(int i = 0; i < size; i++) ptr[i] = data[i];
        return NDArray<T>::from_raw(ptr, shapes);
    }

    template<typename T>
    auto sum(const NDArray<T>& arr)
    {
        return arr.reduce_sum();
    }

    template<typename T>
    T min(const NDArray<T>& arr)
    {
        return arr.reduce_min();
    }

    template<typename T>
    T max(const NDArray<T>& arr)
    {
        return arr.reduce_max();
    }

    template<typename T>
    double mean(const NDArray<T>& arr)
    {
        return arr.reduce_mean();
    }

    template<typename T>
    double median(const NDArray<T>& arr)
    {
        return arr.reduce_median();
    }

    template<typename T>
    double percentile(const NDArray<T>& arr, double p)
    {
        return arr.percentile(p);
    }

    template<typename T>
    double std_(const NDArray<T>& arr)
    {
        return arr.reduce_std();
    }

    template<typename T>
    bool all(const NDArray<T>& arr)
    {
        return arr.reduce_and();
    }

    template<typename T>
    bool any_(const NDArray<T>& arr)
    {
        return arr.reduce_or();
    }

    template<typename T>
    NDArray<T> reshape(const NDArray<T>& arr, const List<int>& new_shape)
    {
        return arr.reshape(new_shape);
    }

    template<typename T>
    NDArray<T> reshape(NDArray<T>&& arr, const List<int>& new_shape)
    {
        return (std::move(arr)).reshape(new_shape);
    }
 
    template<typename T>
    int argmax(const NDArray<T>& arr)
    {
        return arr.argmax();
    }

    template<typename T>
    int argmin(const NDArray<T>& arr)
    {
        return arr.argmin();
    }

    template<typename T>
    NDArray<T> sort(const NDArray<T>& arr, int axis=-1)
    {
        auto res = arr.copy();
        res.sort(axis);
        return res;
    }

    template<typename T>
    auto argsort(const NDArray<T>& arr)
    {
        return arr.argsort();
    }

    template<typename T>
    auto min(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_min(ax);
    }

    template<typename T>
    auto argmin(const NDArray<T>&arr, int ax) 
    {
        return arr.argmin(ax); 
    }

    template<typename T>
    auto max(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_max(ax);
    }

    template<typename T>
    NDArray<int> argmax(const NDArray<T>&arr, int ax) 
    {
        return arr.argmax(ax); 
    }

    template<typename T>
    auto sum(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_sum(ax);
    }

    template<typename T>
    auto mean(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_mean(ax);
    }

    template<typename T>
    auto median(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_median(ax);
    }

    template<typename T>
    auto std_(const NDArray<T>& arr, int ax)
    {
        return arr.reduce_std(ax);
    }


    template<typename T>
    auto all(const NDArray<T>& arr, int ax) 
    {
        return arr.reduce_and(ax);
    }

    template<typename T>
    auto any_(const NDArray<T>&arr, int ax)
    {
        return arr.reduce_or(ax);
    }

    template<typename T>
    auto argsort(const NDArray<T>& arr, int ax=-1)
    {
        return arr.argsort(ax);
    }

    template<typename T>
    auto abs(const NDArray<T>& arr)
    {
        return arr.abs();
    }

    template<typename T>
    auto round(const NDArray<T>& arr)
    {
        return arr.round();
    }

    template<typename T>
    auto ceil(const NDArray<T>& arr)
    {
        return arr.ceil();
    }

    template<typename T>
    auto floor(const NDArray<T>& arr)
    {
        return arr.floor();
    }

    template<typename T>
    auto square(const NDArray<T>& arr)
    {
        return arr.square();
    }

    template<typename T>
    auto sqrt(const NDArray<T>& arr)
    {
        return arr.sqrt();
    }

    template<typename T>
    auto exp(const NDArray<T>& arr)
    {
        return arr.exp();
    }

    template<typename T>
    auto log(const NDArray<T>& arr)
    {
        return arr.log();
    }

    template<typename T>
    auto sin(const NDArray<T>& arr)
    {
        return arr.sin();
    }

    template<typename T>
    auto cos(const NDArray<T>& arr)
    {
        return arr.cos();
    }

    template<typename T>
    auto tan(const NDArray<T>& arr)
    {
        return arr.tan();
    }

    template<typename T>
    auto asin(const NDArray<T>& arr)
    {
        return arr.asin();
    }

    template<typename T>
    auto acos(const NDArray<T>& arr)
    {
        return arr.acos();
    }

    template<typename T>
    auto atan(const NDArray<T>& arr)
    {
        return arr.atan();
    }

    template<typename T>
    auto sinh(const NDArray<T>& arr)
    {
        return arr.sinh();
    }

    template<typename T>
    auto cosh(const NDArray<T>& arr)
    {
        return arr.cosh();
    }

    template<typename T>
    auto tanh(const NDArray<T>& arr)
    {
        return arr.tanh();
    }

    template<typename T>
    auto asinh(const NDArray<T>& arr)
    {
        return arr.asinh();
    }

    template<typename T>
    auto acosh(const NDArray<T>& arr)
    {
        return arr.acosh();
    }

    template<typename T>
    auto atanh(const NDArray<T>& arr)
    {
        return arr.atanh();
    }

    template<typename T>
    NDArray<T> transpose(const NDArray<T>& arr, const List<int>& axes)
    {
        return arr.transpose(axes);
    }

    template<typename T>
    NDArray<T> concatenate(const List<NDArray<T>>& lst, int ax=0)
    {
        return NDArray<T>::concatenate(lst, ax);
    }

    template<typename T>
    NDArray<T> stack(const List<NDArray<T>>& lst, int ax=0)
    {
        return NDArray<T>::stack(lst, ax);
    }

    template<typename T>
    NDArray<T> expand_dims(const NDArray<T>& arr, int ax)
    {
        return arr.expand_dims(ax);
    }

    template<typename T>
    auto minimum(const NDArray<T>& arr, const T& val)
    {
        return arr.template minimum(val);
    }

    template<typename T>
    auto maximum(const NDArray<T>& arr, const T& val)
    {
        return arr.template maximum(val);
    }

    template<typename T>
    auto clip(const NDArray<T>& arr, const T& min_val, const T& max_val)
    {
        return arr.template maximum(min_val).template minimum(min_val);
    }

    template<typename T>
    auto minimum(const NDArray<T>& arr, const NDArray<T>& val)
    {
        return arr.template minimum(val);
    }

    template<typename T>
    auto maximum(const NDArray<T>& arr, const NDArray<T>& val)
    {
        return arr.template maximum(val);
    }

    template<typename T>
    NDArray<T> arange(int stop)
    {
        List<T> lst;
        for(int i = 0; i < stop; i++)
        {
            lst.append(i);
        }
        return array(lst, {stop});
    }

    template<typename T>
    NDArray<T> arange(int start, int stop, Integer step=Integer::none())
    {
        List<T> lst;
        int step_;
        if(!step.is_null)
        {
            step_ = step.value;
        }
        else
        {
            step_ = (start < stop)? 1 : -1;
        }

        assert(step_ != 0);

        if(step_ > 0)
        {
            stop = start + ((stop - start - 1) / step_) * step_ + step_;
        }else
        {
            stop = start - ((start - stop - 1) / (-step_)) * (-step_) + step_;
        }

        
        for(int i = start; i != stop; i += step_)
        {
            lst.append(i);
        }
        return array(lst, {lst.size()});
    }
}

#endif
