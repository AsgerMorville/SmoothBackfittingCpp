//
// Created by Asger Morville on 2024/01/08.
//

#ifndef SMOOTH_BACKFITTING_LIBRARY_QUANTILE_H
#define SMOOTH_BACKFITTING_LIBRARY_QUANTILE_H

#include <algorithm>
#include <cmath>
#include <vector>

template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
    return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs)
{
    if (inData.empty())
    {
        return std::vector<T>();
    }

    if (1 == inData.size())
    {
        return std::vector<T>(1, inData[0]);
    }

    std::vector<T> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<T> quantiles;

    size_t size0 = data.size()-1;
    for (size_t i = 0; i < probs.size(); ++i)
    {
        T pp = size0*probs[i];
        size_t Left = std::floor(pp);
        size_t Right = std::ceil(pp);
        T whole, fractional;
        fractional = std::modf(pp, &whole);
        T datLeft = data.at(Left);
        T datRight = data.at(Right);
        T result = (1-fractional)*datLeft + fractional*datRight;

        quantiles.push_back(result);
    }

    return quantiles;
}

#endif //SMOOTH_BACKFITTING_LIBRARY_QUANTILE_H
