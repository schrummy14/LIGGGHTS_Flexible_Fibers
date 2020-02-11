#ifndef POLYNOMIAL_ROOT_FINDER_H
#define POLYNOMIAL_ROOT_FINDER_H 1

#include <stdio.h>
#include <fstream>
#include <cctype>
#include <cmath>
#include <cfloat>
#include <iostream>

#define MAXDEGREE 100
#define MDP1 MAXDEGREE+1

void rpoly_ak1(double op[MDP1], int* Degree, double zeror[MAXDEGREE], double zeroi[MAXDEGREE]);
void Fxshfr_ak1(int L2, int* NZ, double sr, double v, double K[MDP1], int N, double p[MDP1], int NN, double qp[MDP1], double u, double* lzi, double* lzr, double* szi, double* szr);
void QuadSD_ak1(int NN, double u, double v, double p[MDP1], double q[MDP1], double* a, double* b);
int calcSC_ak1(int N, double a, double b, double* a1, double* a3, double* a7, double* c, double* d, double* e, double* f, double* g, double* h, double K[MDP1], double u, double v, double qk[MDP1]);
void nextK_ak1(int N, int tFlag, double a, double b, double a1, double* a3, double* a7, double K[MDP1], double qk[MDP1], double qp[MDP1]);
void newest_ak1(int tFlag, double* uu, double* vv, double a, double a1, double a3, double a7, double b, double c, double d, double f, double g, double h, double u, double v, double K[MDP1], int N, double p[MDP1]);
void QuadIT_ak1(int N, int* NZ, double uu, double vv, double* szr, double* szi, double* lzr, double* lzi, double qp[MDP1], int NN, double* a, double* b, double p[MDP1], double qk[MDP1], double* a1, double* a3, double* a7, double* c, double* d, double* e, double* f, double* g, double* h, double K[MDP1]);
void RealIT_ak1(int* iFlag, int* NZ, double* sss, int N, double p[MDP1], int NN, double qp[MDP1], double* szr, double* szi, double K[MDP1], double qk[MDP1]);
void Quad_ak1(double a, double b1, double c, double* sr, double* si, double* lr, double* li);

#endif
