#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

#include "modulesCpp/utilityFunctions.h"
#include "modulesCpp/VEC2.h"
#include "modulesCpp/CELL.h"
#include "modulesCpp/constants.h"
#include "modulesCpp/vecCellUtilities.h"
#include "modulesCpp/boundaryFunctions.h"

/*
Units
-----
Distance: micrometres
Time    : hours

Todo: combine the different forces into a single loop
*/

int N;  // #agents
int nt; // #time steps
float rt{70.};
float beta{60.};
float w2, h2;
float W, H;

void writecoordinates(CELL **B, std::ofstream &file);
void writecoordinates(CELL M[], std::ofstream &file);
void sort(float *q, int nelt);

int len(CELL **A);
int in(CELL *p, CELL **m);
int argmax(float M[], int K);

float unitrand();
float map(float, float, float, float, float);
float Hv(float);
float noisemag(float);
float findLargestGap(float *q, int nelt);
float dist(CELL a, CELL b);
float angleAC(VEC2 *rl, VEC2 *rc, VEC2 *rn);

VEC2 mean(CELL **M);
VEC2 mean(CELL M[]);

CELL **getcellarray(int N);
CELL **getNeighbors(CELL M[], CELL *cell, float rt);
CELL **setdiff(CELL **A, CELL **B);

float fInteractionMag(float r)
{
    if (r < rt)
    {
        float U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
        float A0 = 8, A1 = 2, A2 = 25, A3 = 26;
        A0 = 40; // temp
        float force = 0;
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

// void setaccn(CELL &A, CELL &B)
// {
//     VEC2 dp = B.p - A.p;
//     float r = dp.mag();
//     float f = fInteractionMag(r);
//     VEC2 F(f * (dp.x / r), f * (dp.y / r));
//     A.a += F;
// }

VEC2 getInteractionForce(CELL A, CELL B)
{
    VEC2 dp = B.p - A.p;
    float r = dp.mag();
    float f = fInteractionMag(r);
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

/*Loops through every cell and again through every cell  */
void looploop(CELL M[])
{
    for (int i = 0; i < N; i++)
    {
        CELL **B = getNeighbors(M, &M[i], rt);
        VEC2 FIi = getInteractionForce(M[i], B);
        M[i].a = FIi;
    }
    // interaction
    // for (int i = 0; i < N; i++)
    // {
    //     M[i].a = VEC2(0., 0.);
    //     for (int j = 0; j < N; j++)
    //     {

    //         if (i != j)
    //         {
    //             setaccn(M[i], M[j]);
    //         }
    //     }
    // }
    // viscek
    // for (int i = 0; i < N; i++)
    // {
    //     VEC2 a;
    //     int ninr = 0;
    //     VEC2 dv;
    //     for (int j = 0; j < N; j++)
    //     {

    //         VEC2 dp = M[j].p - M[i].p;
    //         float r = dp.mag();
    //         if (r < rt)
    //         {
    //             dv = M[j].v - M[i].v;
    //             a += dv;
    //             ninr += 1;
    //         }
    //     }

    //     M[i].a += a * (beta / ninr);
    // }

    // // noise
    // for (int i = 0; i < N; i++)
    // {
    //     float sig0 = 150.;
    //     float sig1 = 300.;
    //     float rho0 = N / (W * H); // reference density (change to rho1)
    //     float rho = 0.;
    //     float sig = 0.;
    //     float tau = 1.39;
    //     // float tau = 0.01;

    //     VEC2 U = VEC2(randf(0, 1), randf(0, 1));
    //     VEC2 xi(sqrt(-2 * log(U.x)) * cos(2 * PI * U.y),
    //             sqrt(-2 * log(U.x)) * sin(2 * PI * U.y));
    //     float std1 = 1; // temp
    //     xi = xi * std1;
    //     float theta = 0.;
    //     for (int j = 0; j < N; j++)
    //     {

    //         VEC2 dp = M[j].p - M[i].p;
    //         float r = dp.mag();
    //         if (r < rt)
    //         {
    //             rho += 1;
    //         }
    //     }
    //     rho /= PI * rt * rt;
    //     sig = sig0 + (sig1 - sig0) * (1 - rho / rho0);
    //     // euler maruyama integration
    //     M[i].eta = M[i].eta - (M[i].eta) * (dt / tau);
    //     M[i].eta = M[i].eta + (xi) * (sqrt(dt) / tau);

    //     M[i].eta = M[i].eta / M[i].eta.mag(); // normalize eta
    //     M[i].a += M[i].eta * sig;
    // }

    for (int i = 0; i < N; i++)
    {
        M[i].update();
    }
}

int main(int argc, char *argv[])
{
    N = std::stoi(argv[1]);
    nt = std::stoi(argv[2]);
    w2 = W / 2;
    h2 = H / 2;
    rt = std::stof(argv[3]);

    std::fstream pointinshape;
    pointinshape.open("intermediateResults/pointsInPolygon.csv", std::ios_base::in);

    CELL *M = (CELL *)malloc(sizeof(CELL) * N);
    float x, y;
    for (int i = 0; i < N; i++)
    {
        pointinshape >> x >> y;
        M[i] = CELL(x, y);
    }

    std::ofstream posfile;
    posfile.open("intermediateResults/positionData.csv");

    // std::ofstream boundposfile;
    // boundposfile.open("intermediateResults/boundaryData.csv");
    // CELL **B;

    for (int t = 0; t < nt; t++)
    {
        looploop(M);
        writecoordinates(M, posfile);
        // B = findBorderCellsByEdgeScan(M, 200);
        // writecoordinates(B, boundposfile);
    }

    std::cout << "simulation complete\n";
}