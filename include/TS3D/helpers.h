
#ifndef common_h
#define common_h

#include <iostream>
#include <limits>
#include <set>

#include <Eigen/Sparse>

#include "timer.h"

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double, RowMajor> SparseMatrixd;
typedef SparseMatrix<int, RowMajor> SparseMatrixi;

void divideRowByGCD(SparseMatrixi& C, int i);
void divideRowByGCD(SparseMatrixi& C);

void divideRow(SparseMatrixd& C, int i);

void IREF(SparseMatrixi& C, VectorXd& b, VectorXi& indexC, VectorXi& indexR, bool gcd_division, double& total_time);
void IRREF(SparseMatrixi& C, VectorXd& b, VectorXi indexC, bool gcd_division, double& total_time);
VectorXd
evaluate(SparseMatrixi C, VectorXd X_bar, VectorXd b, double M, double uv_max, VectorXi indexC, VectorXi indexR);

double safeDotProd(vector<pair<int, double>> S, double M);
double makeDiv(double x, vector<int> D, double& f_max);
double fixedPrec(double x, double& f_max);
double fixedPrec(VectorXd& X_bar, int k, double& f_max);

int gcd(int a, int b);
int lcm(int a, int b);
int lcmArr(vector<int> D);

void statistic(SparseMatrixi C);
void statistic(SparseMatrixi C, VectorXi indexC);

int divideInt(int dividend, int divisor, int s = 0);
int addInt(int a, int b);
int multiplyInt(int a, int b);

#endif /* common_h */
