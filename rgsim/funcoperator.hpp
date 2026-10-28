#pragma once
#include "types.hpp"

namespace rgsim {

inline func constant(double v)
{
    return [v](const std::array<dual,N>&) { dual res = v; return res; };
}

inline func operator+(func a, func b)
{
    return [a, b](const std::array<dual,N>& x)
    { dual p = a(x); dual q = b(x); dual res = p + q; return res; };
}
inline func operator-(func a, func b)
{
    return [a, b](const std::array<dual,N>& x)
    { dual p = a(x); dual q = b(x); dual res = p - q; return res; };
}
inline func operator*(func a, func b)
{
    return [a, b](const std::array<dual,N>& x)
    { dual p = a(x); dual q = b(x); dual res = p * q; return res; };
}
inline func operator/(func a, func b)
{
    return [a, b](const std::array<dual,N>& x)
    { dual p = a(x); dual q = b(x); dual res = p / q; return res; };
}

} // namespace rgsim