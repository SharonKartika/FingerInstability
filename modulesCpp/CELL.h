#ifndef CELL_H
#define CELL_H
#include "VEC2.h"
#include "constants.h"
class CELL
{
public:
    VEC2 p, v, a, eta;
    CELL()
    {
        float theta = randf(0, 2 * PI);
        float r = w2 * sqrt(randf(0, 1));
        p = VEC2(r * cos(theta), r * sin(theta));
        // p = VEC2(randf(-w2, w2), randf(-h2, h2));
        v = VEC2(randf(-10, 10), randf(-10, 10));
        a = VEC2(0., 0.);
        theta = randf(0, 2 * PI);
        eta = VEC2(cos(theta), sin(theta));
    }
    CELL(float x, float y)
    {
        float theta = randf(0, 2 * PI);
        p = VEC2(x, y);
        v = VEC2(randf(-10, 10), randf(-10, 10));
        a = VEC2(0., 0.);
        theta = randf(0, 2 * PI);
        eta = VEC2(cos(theta), sin(theta));
    }
    void update()
    {
        v = v + a * dt;
        p = p + v * dt;
    }
};
#endif
