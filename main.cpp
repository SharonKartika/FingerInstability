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
#include "modulesCpp/interactionForce.h"

/*
Units
-----
Distance: micrometres
Time    : hours

Todo: combine the different forces into a single loop
*/

int N;  // #agents
int nt; // #time steps
double rt{70.};
double beta{60.};
double w2, h2;
double W, H;

void writecoordinates(CELL **B, std::ofstream &file);
void writecoordinates(CELL M[], std::ofstream &file);
void sort(double *q, int nelt);

int len(CELL **A);
int in(CELL *p, CELL **m);
int argmax(double M[], int K);

double unitrand();
double map(double, double, double, double, double);
double Hv(double);
double noisemag(double);
double findLargestGap(double *q, int nelt);
double dist(CELL a, CELL b);
double angleAC(VEC2 *rl, VEC2 *rc, VEC2 *rn);

VEC2 mean(CELL **M);
VEC2 mean(CELL M[]);

CELL **getcellarray(int N);
CELL **getNeighbors(CELL M[], CELL *cell, double rt);
CELL **setdiff(CELL **A, CELL **B);

// interaction force functions
double fInteractionMag(double r);
VEC2 getInteractionForce(CELL A, CELL B);
VEC2 getInteractionForce(CELL A, CELL **B);

VEC2 getVicsekForce(CELL A, CELL **B)
{
    VEC2 FVc(0.0, 0.0);
    int ninr = 0;
    VEC2 dv(0.0, 0.0);
    while (*B)
    {
        dv = ((*B)->v - A.v);
        FVc += dv;
        ninr++;
        B++;
    }
    if (ninr == 0)
        return VEC2(0.0, 0.0);
    return FVc / double(ninr);
}

VEC2 unitnoise()
{
    VEC2 U = VEC2(randf(0, 1), randf(0, 1));
    VEC2 xi(sqrt(-2 * log(U.x)) * cos(2 * PI * U.y),
            sqrt(-2 * log(U.x)) * sin(2 * PI * U.y));
    return xi;
}

double getLocalCellDensity(CELL **B)
{
    double rho = len(B) + 1; // including the current cell
    return rho / (PI * rt * rt);
}

VEC2 getNoiseForce(CELL &A, CELL **B)
{
    double sig0 = 150.,
           sig1 = 300.,
           rho0 = N / (W * H),
           rho = 0.,
           sig = 0.,
           tau = 1.39;
    VEC2 xi = unitnoise();
    rho = getLocalCellDensity(B);
    sig = sig0 + (sig1 - sig0) * (1 - rho / rho0);

    // integrate and normalize
    A.eta -= A.eta * dt / tau;
    A.eta += xi * sqrt(dt) / tau;
    A.eta = A.eta / A.eta.mag();

    return A.eta * sig;
}

/*Loops through every cell and again through every cell  */
void looploop(CELL M[])
{
    for (int i = 0; i < N; i++)
    {
        CELL **B = getNeighbors(M, &M[i], rt);
        VEC2 FIn, FVc, FNo;

        FIn = getInteractionForce(M[i], B);
        FVc = getVicsekForce(M[i], B);
        FNo = getNoiseForce(M[i], B);

        M[i].a = FIn +
                 FVc * beta +
                 FNo;
    }
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
    double x, y;
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