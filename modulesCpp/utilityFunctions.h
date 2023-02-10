#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>
#include "constants.h"

double unitrand()
/*Returns a double chosen randomly from [0,1]  */
{
    return double(rand()) / INT_MAX;
}

double map(double ri, double x1, double x2, double y1, double y2)
/*Maps one interval to another */
{
    double runit = (ri - x1) / (x2 - x1);
    double rf = runit * (y2 - y1) + y1;
    return rf;
}

double randf(double a = 0., double b = 1.)
/*Returns a double chosen randomly
Defaults to [0,1]*/
{
    double r = double(rand()) / INT_MAX;
    return map(r, 0, 1, a, b);
}


double Hv(double r)
{ // heaviside function
    return r > 0;
}

double noisemag(double rho)
{ // noise magnitude
    double s0 = 150;
    double s1 = 300;
    double rho0 = 2.2e-3;
    return s0 + (s1-s0)*(1-(rho/rho0));
}


/* Takes in an array and its lengths.
Finds index of largest element.*/
int argmax(double M[], int K)
{
    int mi = 0;
    for (int i = 1; i < K; i++)
    {
        if (M[mi] <= M[i])
        {
            mi = i;
        }
    }
    return mi;
}

double findLargestGap(double *q, int nelt)
{
    double lg = 0.;
    double t;
    for (int i = 0; i < nelt - 1; i++)
    {
        t = q[i + 1] - q[i];
        if (t > lg)
        {
            lg = t;
        }
    }
    t = 2 * PI - q[nelt - 1] + q[0];
    if (t > lg)
    {
        lg = t;
    }
    return lg;
}

void sort(double *q, int nelt)
{
    double key;
    int j;
    for (int i = 1; i < nelt; i++)
    {
        key = q[i];
        j = i - 1;

        while (j >= 0 && q[j] > key)
        {
            q[j + 1] = q[j];
            j = j - 1;
        }
        q[j + 1] = key;
    }
}

#endif