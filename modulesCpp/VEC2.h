#ifndef VEC2_H
#define VEC2_H
#include "utilityFunctions.h"

class VEC2
{
public:
    double x, y;
    VEC2(double X, double Y)
    {
        x = X;
        y = Y;
    }
    VEC2()
    {
        x = 0.;
        y = 0.;
    }
    VEC2 operator+(VEC2 const &obj)
    {
        VEC2 temp;
        temp.x = x + obj.x;
        temp.y = y + obj.y;
        return temp;
    }
    VEC2 operator-(VEC2 const &obj)
    {
        VEC2 temp;
        temp.x = x - obj.x;
        temp.y = y - obj.y;
        return temp;
    }
    VEC2 operator*(double const &a)
    {
        VEC2 temp;
        temp.x = a * x;
        temp.y = a * y;
        return temp;
    }
    void operator+=(VEC2 const &obj)
    {
        x += obj.x;
        y += obj.y;
    }
    void operator-=(VEC2 const &obj)
    {
        x -= obj.x;
        y -= obj.y;
    }
    VEC2 operator/(double const &a)
    {
        VEC2 temp;
        temp.x = x / a;
        temp.y = y / a;
        return temp;
    }
    double mag()
    {
        return sqrt(pow(x, 2) + pow(y, 2));
    }
    VEC2 unit(){
        return (*this)/this->mag();
    }
};

#endif
