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
double W = 1200, H = 1200;

std::fstream lbc; // lines between 2 closest

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

/*returns array of indices of k smallest elements
Destroys the input array in the process*/
int *argkmin(double q[], int n, int k)
{
    int *kq = (int *)malloc(sizeof(int) * k);
    for (int i = 0; i < k; i++)
    {
        int j = argmin(q, n);
        q[j] = std::numeric_limits<double>::max();
        kq[i] = j;
    }
    return kq;
}

CELL **kClosest(CELL *T, CELL **B, int k)
{
    CELL **KN = getCellPointerArray(k);
    int n = len(B);
    double *dists = (double *)malloc(sizeof(double) * n);
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
    VEC2 dp = (B->p - A->p);
    double r = dp.mag();
    VEC2 u = dp.unit();
    double f = 1 / (r * r);
    return u * f;
}

void writeKBoundInteraction(CELL *T, CELL *A)
{
    lbc << T->p.x << ","
        << T->p.y << ","
        << A->p.x << ","
        << A->p.y << ",";
}

VEC2 testAttractiveForce(CELL *A, CELL *B)
{
    writeKBoundInteraction(A, B);

    VEC2 dp = (B->p - A->p);
    double r = dp.mag();
    VEC2 unit = dp.unit();

    //Force 1
    // double k = 100;
    // VEC2 F;
    // double rcutoff = 300.0;
    // if (r > rcutoff)
    //     return VEC2(0, 0);
    // else
    //     return unit * r * k;
    
    //Force 2
    // double mag = (r > rcutoff) ? 0.0 : -(r - rcutoff);
    // F = unit * mag * k;
    // return F;

    //Force 3
    return unit * -(1/pow(r,3)-10)*100;
}

// interaction force between boundary cells
VEC2 boundaryInteractionForce(CELL *A, CELL *B)
{
}

// calculate force on T due to A and B
VEC2 getActinForce(CELL *T, CELL *A, CELL *B)
{
    VEC2 FAc = testAttractiveForce(T, A) +
               testAttractiveForce(T, B);
    return FAc;
}

// void writeKBoundInteraction(CELL *T, CELL *A, CELL *B)
// {
//     lbc << T->p.x << ","
//         << T->p.y << ","
//         << A->p.x << ","
//         << A->p.y << ","
//         << T->p.x << ","
//         << T->p.y << ","
//         << B->p.x << ","
//         << B->p.y << ",";
// }



VEC2 getActinForce(CELL *T, CELL **B)
{
    CELL **KN = kClosest(T, B, 2);
    VEC2 FAc = getActinForce(T, *(KN + 0), *(KN + 1));
    // writeKBoundInteraction(T, *(KN + 0), *(KN + 1));
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
        VEC2 FIn;
        VEC2 FVc;
        VEC2 FNo;

        FIn = getInteractionForce(M[i], B);
        FVc = getVicsekForce(M[i], B);
        FNo = getNoiseForce(M[i], B);

        M[i].a =// VEC2(0, 0);
            FIn +
            FVc * beta +
            FNo;
    }
    CELL **B = findBorderCellsByFOV(M, 400, PI / 2);
    for (int i = 0; i < len(B); i++)
    {

        //attract the closest two cells
        // VEC2 FAc; // F bodiesactin
        // FAc = getActinForce(*(B + i), B);
        // (*(B + i))->a += FAc;
        
        //attract within a neighborhood, on the boundary
        CELL **D = getNeighbors(B, *(B + i), 300);
        for (int j = 0; j < len(D); j++)
            (*(B + i))->a += testAttractiveForce(*(B + i), *(D + j));
        
        //add noise to boundarycells, due to ALL neighboring cells
        // CELL **D = getNeighbors(M, *(B + i), rt);
        // (*(B + i))->a += getNoiseForce(*(*(B + i)), D);
    }
    lbc << std::endl;
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

    lbc.open("intermediateResults/linesBetweenCells.csv", std::fstream::out);
    lbc.close();
    lbc.open("intermediateResults/linesBetweenCells.csv", std::fstream::app);
    std::fstream pointinshape;
    pointinshape.open("intermediateResults/pointsInPolygon.csv", std::ios_base::in);

    CELL M[N];
    // CELL *M = (CELL *)malloc(sizeof(CELL) * N);
    // double x, y;
    // for (int i = 0; i < N; i++)
    // {
    //     pointinshape >> x >> y;
    // M[i] = CELL(x, y);
    // M[i] = CELL();
    // }

    std::ofstream posfile;
    posfile.open("intermediateResults/positionData.csv");

    std::ofstream boundposfile;
    boundposfile.open("intermediateResults/boundaryData.csv");

    for (int t = 0; t < nt; t++)
    {
        looploop(M, posfile, boundposfile);
    }
    lbc.close();
    std::cout << "simulation complete\n";
}