#ifndef BOUNDARYFUNCTIONS_H
#define BOUNDARYFUNCTIONS_H
#include "utilityFunctions.h"
#include "VEC2.h"
#include "CELL.h"
#include "constants.h"
#include "vecCellUtilities.h"

float angleAC(VEC2 *rl, VEC2 *rc, VEC2 *rn)
{
    VEC2 A = *rc - *rl;
    VEC2 B = *rn - *rc;
    float k = B.y * A.x - A.y * B.x;
    float costheta = (A.x * B.x + A.y * B.y) /
                     (A.mag() * B.mag());
    float theta = acos(costheta);
    if (k > 0)
        return PI - theta;
    else
        return 2 * PI - theta;
}


/*Writes to a file the cells that lie on the border
Takes the farthest cell.
Finds the next cell on the boundary incrementally*/
CELL **findBorderCellsByEdgeScan(CELL M[], float rt)
{
    VEC2 rcm = mean(M);
    float dists[N];
    for (int i = 0; i < N; i++)
    {
        dists[i] = (rcm - M[i].p).mag();
    }
    int i = argmax(dists, N);
    CELL *rft = &M[i];
    CELL *rc = rft;

    CELL **boundcells = getcellarray(N);
    while (true)
    {
        CELL **rns = getNeighbors(M, rc, rt);
        VEC2 rl = mean(rns);
        CELL **rnsnb = setdiff(rns, boundcells);
        int lenrnsnb = len(rnsnb);
        float *nscores = (float *)malloc(lenrnsnb * sizeof(*nscores));
        for (int i = 0; i < lenrnsnb; i++)
        {
            *(nscores + i) = angleAC(&rl, &(rc->p), &(*(rnsnb + i))->p);
            if (*(nscores + i) >= PI)
                (*(nscores + i) = 0.);
        }
        // debug segfault
        std::cout << len(rns) << "\t" << len(rnsnb) << "\t" << len(boundcells) << std::endl;
        // end debug segfault

        int rnindex = argmax(nscores, lenrnsnb);
        CELL *rn = rnsnb[rnindex];
        rc = rn;
        boundcells[len(boundcells)] = rc;

        if (rc == rft)
        {
            std::cout << "DONE ";
            break;
        }
    }
    return boundcells;
}

bool isQuadrantEmpty(int *a)
{
    for (int i = 0; i < 4; i++)
    {
        if (a[i] == 0)
        {
            return true;
        }
    }
    return false;
}

int *getQuadrantCount(CELL **NC, CELL *C)
{
    static int qc[4];
    for (int i = 0; i < 4; i++)
    {
        qc[i] = 0;
    }
    int Yv, Xv;
    int counttemp = 0;
    while (*NC)
    {
        if ((*NC)->p.x > C->p.x)
            Xv = 1;
        else
            Xv = 0;
        if ((*NC)->p.y > C->p.y)
            Yv = 1;
        else
            Yv = 0;

        if (Xv == 1 && Yv == 1)
            qc[0]++;
        else if (Xv == 1 && Yv == 0)
            qc[1]++;
        else if (Xv == 0 && Yv == 0)
            qc[2]++;
        else if (Xv == 0 && Yv == 1)
            qc[3]++;
        else
            std::cout << "Some kind of error" << std::endl;
        NC++;
    }
    return qc;
}

bool isOnBoundary(int *a)
{
    // Difference between number of cells in one region
    //  greater than 5
    //  than the average number of cells in other 3
    for (int i = 0; i < 4; i++)
    {
        float avg = 0;
        for (int j = 0; j < 4; j++)
        {
            if (i != j)
            {
                avg += a[j];
            }
        }
        avg /= 3;
        if (abs(a[i] - avg) > 5)
        {
            return true;
        }
    }
    // Difference between total particles in 2 regions
    //  greater than 8
    //  than total particles in the other 2
    for (int i = 0; i < 4; i++)

    {
        int c1, c2;
        c1 = a[i % 4] + a[(i + 1) % 4];
        c2 = a[(i + 2) % 4] + a[(i + 3) % 4];
        if (abs(c1 - c2) > 8)
        {
            return true;
        }
    }
    // opposite quadrants, count > 8?
    int v1 = a[0] + a[2];
    int v2 = a[1] + a[3];
    if (abs(v1 - v2) > 8)
    {
        return true;
    }
    return false;
}

// -[ ] Optimise using short circuiting
bool isOnBoundaryFOV(CELL **rns, CELL *C, float f)
{
    int nelt = len(rns);
    float *q = (float *)malloc(sizeof(float) * nelt);
    for (int i = 0; i < nelt; i++)
    {
        VEC2 s = (*(rns + i))->p - C->p;
        *(q + i) = atan2(s.y, s.x);
    }
    sort(q, nelt);
    float lg = findLargestGap(q, nelt);
    if (lg > f)
    {
        return true;
    }
    return false;
}

CELL **findBorderCellsByFOV(CELL M[], float rt, float f)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        CELL **rns = getNeighbors(M, &M[i], rt);
        if (isOnBoundaryFOV(rns, &M[i], f))
        {
            *p = &M[i];
            p++;
        }
    }
    std::cout << "found by FOV" << std::endl;
    return boundcells;
}

CELL **findBorderCellsByVecSum(CELL M[], float rt, float trmag)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        VEC2 vecsum = VEC2(0, 0);
        CELL **rns = getNeighbors(M, &M[i], rt);
        CELL **t = rns;
        while (*t)
        {
            vecsum += ((*t)->p - M[i].p);
            t++;
        }
        vecsum = vecsum / len(rns);
        if (vecsum.mag() > trmag)
        {
            std::cout << vecsum.mag() << std::endl;
            *p = &M[i];
            p++;
        }
    }
    std::cout << "found by Vecsum" << std::endl;
    return boundcells;
}

CELL **findBorderCellsByLevine(CELL M[], float rt)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        CELL **rns = getNeighbors(M, &M[i], rt);
        int *qc = getQuadrantCount(rns, &M[i]);
        if (isOnBoundary(qc))
        {
            *p = &M[i];
            p++;
        }
    }
    std::cout << "Found by levine method" << std::endl;
    return boundcells;
}

CELL **findBorderCellsByQuadrantEmpty(CELL M[], float rt)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        CELL **rns = getNeighbors(M, &M[i], rt);
        int *qc = getQuadrantCount(rns, &M[i]);
        if (isQuadrantEmpty(qc))
        {
            *p = &M[i];
            p++;
        }
    }
    std::cout << "Found by quadrantempty" << std::endl;
    return boundcells;
}

#endif