#ifndef INTERACTIONFORCE_H
#define INTERACTIONFORCE_H

#include <math.h>
#include "constants.h"
#include "VEC2.h"
#include "CELL.h"

double fInteractionMag(double r)
{
    if (r < rt)
    {
        double U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
        double A0 = 8, A1 = 2, A2 = 25, A3 = 26;
        double force = 0;
        force += U0 * r * exp(-(pow((r / A0), 2)));
        force += U2 * exp(-r / A2);
        force -= U3 * pow(r - A3, 2) * Hv(r - A3);
        force += U1 * (r - A1) * Hv(r - A1);
        return -force;
    }
    else
    {
        return 0.;
    }
}

VEC2 getInteractionForce(CELL A, CELL B)
{
    VEC2 dp = B.p - A.p;
    double r = dp.mag();
    double f = fInteractionMag(r);
    VEC2 F(f * (dp.x / r), f * (dp.y / r));
    return F;
}

VEC2 getInteractionForce(CELL A, CELL **B)
{
    VEC2 FIi(0.0, 0.0); // interaction force on i
    while (*B)
    {
        FIi += getInteractionForce(A, **B);
        B++;
    }
    return FIi;
}

#endif