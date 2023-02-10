#ifndef VECCELLUTILITIES_H
#define VECCELLUTILITIES_H

#include "utilityFunctions.h"
#include "VEC2.h"
#include "CELL.h"
#include "constants.h"

/*Update this. There are much more elegant ways*/
void cout(VEC2 *p)
{
    std::cout << p->x << "  " << p->y << "  ";
}

int len(CELL **A)
{
    int i = 0;
    while (*A)
    {
        A++;
        i++;
    }
    return i;
}

void writecoordinates(CELL M[], std::ofstream &file)
{ /*
  Writes the coordinates of all the cells to file
  */
    for (int i = 0; i < N - 1; i++)
        file << M[i].p.x << ',' << M[i].p.y << ',';
    file << M[N - 1].p.x << ',' << M[N - 1].p.y << '\n';
}

void writecoordinates(CELL **B, std::ofstream &file)
{
    int nb = len(B);
    for (int i = 0; i < nb - 1; i++)
        file << (*(B + i))->p.x << ',' << (*(B + i))->p.y << ',';
    file << (*(B + nb - 1))->p.x << ',' << (*(B + nb - 1))->p.y << '\n';
}

VEC2 mean(CELL **M)
{
    VEC2 m(0, 0);
    CELL **p = M;
    float i = 0.;

    while (*p != NULL)
    {
        m += (*p)->p;
        i += 1.;
        p++;
    }
    return m / i;
}

VEC2 mean(CELL M[])
{
    VEC2 m(0, 0);
    for (int i = 0; i < N; i++)
    {
        m += M[i].p;
    }
    return m / N;
}

/*  Allocate N pointers.
    Return pointer to pointer array
    */
CELL **getcellarray(int N)
{
    CELL **p = (CELL **)malloc(N * sizeof(p));
    for (int i = 0; i < N; i++)
    {
        p[i] = NULL;
    }
    return p;
}


float dist(CELL a, CELL b)
{
    return (a.p - b.p).mag();
}


CELL **getNeighbors(CELL M[], CELL *cell, float rt)
{
    CELL **rns = getcellarray(N);

    CELL **p = rns;
    for (int i = 0; i < N; i++)
    {
        if (dist(*cell, M[i]) < rt)
        {
            if (cell != &M[i])
            {
                *p = &M[i];
                p++;
            }
        }
    }
    return rns;
}

/*Checks if A is in B*/
int in(CELL *p, CELL **m)
{
    while (*m)
    {
        if (p == *m)
        {
            return 1;
        }
        m++;
    }
    return 0;
}

/*  Get set difference A\B
    Elements in A, not in B.
    */
CELL **setdiff(CELL **A, CELL **B)
{
    CELL **p = getcellarray(N / 5);
    CELL **c = p;
    while (*A)
    {
        if (!in(*A, B))
        {
            *c = *A;
            c++;
        }
        A++;
    }
    return p;
}

#endif