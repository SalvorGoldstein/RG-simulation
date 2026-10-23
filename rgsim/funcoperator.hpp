#pragma once
#include "types.hpp"

namespace rgsim {

inline func constant(double v)
{
    return [v](dual, dual, dual, dual) { dual res = v; return res; };
}

inline func operator+(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x + y; return res; };
}
inline func operator-(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x - y; return res; };
}
inline func operator*(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x * y; return res; };
}
inline func operator/(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x / y; return res; };
}

} // namespace rgsim