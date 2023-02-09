#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>
#include "constants.h"

float unitrand()
/*Returns a float chosen randomly from [0,1]  */
{
    return float(rand()) / INT_MAX;
}

float map(float ri, float x1, float x2, float y1, float y2)
/*Maps one interval to another */
{
    float runit = (ri - x1) / (x2 - x1);
    float rf = runit * (y2 - y1) + y1;
    return rf;
}

float randf(float a = 0., float b = 1.)
/*Returns a float chosen randomly
Defaults to [0,1]*/
{
    float r = float(rand()) / INT_MAX;
    return map(r, 0, 1, a, b);
}


float Hv(float r)
{ // heaviside function
    return r > 0;
}

float noisemag(float rho)
{ // noise magnitude
    float s0 = 150;
    float s1 = 300;
    float rho0 = 2.2e-3;
    return s0 + (s1-s0)*(1-(rho/rho0));
}


/* Takes in an array and its lengths.
Finds index of largest element.*/
int argmax(float M[], int K)
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

float findLargestGap(float *q, int nelt)
{
    float lg = 0.;
    float t;
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

void sort(float *q, int nelt)
{
    float key;
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