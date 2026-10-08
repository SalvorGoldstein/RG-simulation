#include "funcoperator.h"

using namespace std;

func constant(double v)
{
    return [v](dual t, dual r, dual th, dual ph) { dual res = v; return res; };
}

func operator+(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x + y; return res; };
}
func operator-(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x - y; return res; };
}
func operator*(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x * y; return res; };
}
func operator/(func a, func b)
{
    return [a, b](dual t, dual r, dual th, dual ph)
    { dual x = a(t,r,th,ph); dual y = b(t,r,th,ph); dual res = x / y; return res; };
}