#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>
#include <limits>

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

CELL **getCellPointerArray(int N);
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

int argmin(double q[], int n)
{
    int key = 0;
    for (int i = 1; i < n; i++)
        if (q[i] < q[key])
            key = i;
    return key;
}

void swap(auto &a, auto &b)
{
    auto t = a;
    a = b;
    b = t;
}

int *argkmin(double q[], int n, int k)
{
    int *kq = (int *)malloc(sizeof(int) * k);
    for (int i = 0; i < k; i++)
    {
        int j = argmin(q + i, n - i);
        swap(q[i], q[j]);
        kq[i] = j;
    }
    return kq;
}

CELL **kClosest(CELL *T, CELL **B, int k)
{
    CELL **KN = getCellPointerArray(k);
    int n = len(B);
    double *dists = (double *)malloc(sizeof(double) * len(B));
    for (int i = 0; i < n; i++)
    {
        dists[i] = ((*(B + i))->p - T->p).mag();
        if (*(B + i) == T) // if equal, set infinite distance
            dists[i] = std::numeric_limits<double>::max();
    }
    int *km = argkmin(dists, n, k);
    for (int i = 0; i < k; i++)
    {
        *(KN + i) = *(B + km[i]);
    }
    return KN;
}

VEC2 getGravityForce(CELL *A, CELL *B)
{
    VEC2 dp = (A->p - B->p);
    double r = dp.mag();
    double fmag = 1 / (r * r);
    VEC2 FAc = VEC2(fmag * dp.x / r,
                    fmag * dp.y / r);
    return FAc;
}
// calculate force on T due to A and B
VEC2 getActinForce(CELL *T, CELL *A, CELL *B)
{
    // VEC2 FAc = getInteractionForce(*T, *A) +
    //            getInteractionForce(*T, *B);
    // return FAc * 10;
    VEC2 FAc = getGravityForce(T, A) + getGravityForce(T, B);
    const double scale = 1E6;
    // std::cout << FAc.x << " " << FAc.y << std::endl;
    // exit(0);
    return FAc * scale;
}

VEC2 getActinForce(CELL *T, CELL **B)
{
    CELL **KN = kClosest(T, B, 2);
    VEC2 FAc = getActinForce(T, *(KN + 0), *(KN + 1));
    return FAc;
}

/*Loops through every cell and again through every cell  */
void looploop(CELL M[],
              std::ofstream &posfile,
              std::ofstream &boundposfile)
{
    for (int i = 0; i < N; i++)
    {
        CELL **B = getNeighbors(M, &M[i], rt);
        VEC2 FIn, FVc, FNo;

        FIn = getInteractionForce(M[i], B);
        FVc = getVicsekForce(M[i], B);
        FNo = getNoiseForce(M[i], B);

        M[i].a = VEC2(0, 0);
        //  FIn +
        //  FVc * beta +
        //  FNo;
    }
    CELL **B = findBorderCellsByFOV(M, 400, 3 * PI / 4);
    for (int i = 0; i < len(B); i++)
    {
        VEC2 FAc; // F actin
        FAc = getActinForce(*B, B);
        (*(B + i))->a += FAc;
    }
    writecoordinates(M, posfile);
    writecoordinates(B, boundposfile);
    // boundposfile.close();
    // exit(1);
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

    std::ofstream boundposfile;
    boundposfile.open("intermediateResults/boundaryData.csv");

    

    // for (int t = 0; t < nt; t++)
    // {
        // looploop(M, posfile, boundposfile);
    // }

    std::cout << "simulation complete\n";
}