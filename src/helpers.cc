#include "TS3D/helpers.h"

#include <math.h>
#include <stack>

void divideRowByGCD(SparseMatrixi& C)
{
    for (int i = 0; i < C.outerSize(); i++)
        divideRowByGCD(C, i);
}

void divideRowByGCD(SparseMatrixi& C, int i)
{
    int min_val = -1;
    for (SparseMatrixi::InnerIterator it(C, i); it; ++it)
    {
        if (it.value() != 0)
        {
            int v = it.value();

            if (v < 0)
                v = v * (-1);
            if (min_val < 0)
                min_val = v;
            else
                min_val = gcd(min_val, v);
        }
    }

    if (min_val > 1)
    {
        for (SparseMatrixi::InnerIterator it(C, i); it; ++it)
        {
            if (it.value() != 0)
                C.coeffRef(i, it.index()) = it.value() / min_val;
        }
    }
}

void divideRow(SparseMatrixd& C, int i)
{
    double min_val = -1.0;
    for (SparseMatrixd::InnerIterator it(C, i); it; ++it)
    {
        if (it.value() != 0)
        {
            int v = it.value();

            if (v < 0)
                v = v * (-1.0);
            if (min_val < 0 || min_val < v)
                min_val = v;
        }
    }

    if (min_val > 1.0)
    {
        for (SparseMatrixd::InnerIterator it(C, i); it; ++it)
        {
            if (it.value() != 0)
                C.coeffRef(i, it.index()) = it.value() / min_val;
        }
    }
}

/* Algorithm 1 IREF by fraction-free Gaussian Elimination */
void IREF(SparseMatrixi& C, VectorXd& b, VectorXi& indexC, VectorXi& indexR, bool gcd_division, double& total_time)
{
    auto start = std::chrono::high_resolution_clock::now();

    indexC = VectorXi::Constant(C.outerSize(), -1); // column index of the pivot
    indexR = VectorXi::Constant(C.innerSize(), -1); // row index of the pivot

    int dividend = 0;
    int o = 1;
    int r = 0;

    for (int k = 0; k < C.outerSize(); k++)
    {
        while (k + r < C.cols())
        {
            if (C.coeff(k, k + r) != 0)
            {
                indexC(k) = k + r;
                indexR(k + r) = k;
                break;
            }

            bool undefined = true; // column is all zeroes
            for (int l = k + 1; l < C.outerSize(); l++)
            {
                if (C.coeff(l, k + r) != 0)
                {
                    // find pivot in column & swap rows
                    SparseVector<int> tmp_k = C.row(k).transpose();
                    SparseVector<int> tmp_l = C.row(l).transpose();
                    C.row(k) = tmp_l.transpose();
                    C.row(l) = tmp_k.transpose();
                    undefined = false;
                    break;
                }
            }

            if (undefined)
                r++;
        }

        if (k + r >= C.cols())
            break;

        if (!gcd_division)
            divideRowByGCD(C, k);

        // elimination(row i by j)
        int Ckkr = C.coeff(k, k + r);

        for (int i = k + 1; i < C.outerSize(); i++)
        {
            SparseMatrixi::InnerIterator it_i(C, i);
            int Cikr = C.coeff(i, k + r);

            if (Cikr == 0)
            {
                for (; it_i; ++it_i)
                {
                    if (it_i.index() <= k + r && it_i.value() != 0)
                        it_i.valueRef() = 0;

                    if (!gcd_division)
                    {
                        if (it_i.index() > k + r && it_i.value() != 0)
                        {
                            dividend = multiplyInt(Ckkr, it_i.value());
                            it_i.valueRef() = divideInt(dividend, o);
                        }
                    }
                }
            }
            else
            {
                int local_gcd = -1;
                vector<pair<int, int>> upd;
                SparseMatrixi::InnerIterator it_k(C, k);

                while (it_i || it_k)
                {
                    while ((it_i) && it_i.value() == 0)
                        ++it_i;
                    while ((it_k) && it_k.value() == 0)
                        ++it_k;

                    if (!(it_i || it_k))
                        break;

                    int ix = (it_i) ? it_i.index() : C.cols();
                    int kx = (it_k) ? it_k.index() : C.cols();
                    int j = min(ix, kx);

                    if (j > k + r)
                    {
                        int Cij = (ix == j) ? it_i.value() : 0;
                        int Ckj = (kx == j) ? it_k.value() : 0;

                        dividend = addInt(multiplyInt(Ckkr, Cij), -multiplyInt(Cikr, Ckj));
                        int v = divideInt(dividend, o);
                        upd.push_back(pair<int, int>(j, v));
                        if (gcd_division && 0 != v)
                        {
                            if (v < 0)
                                v *= -1;
                            if (local_gcd < 0)
                                local_gcd = v;
                            else
                                local_gcd = gcd(local_gcd, v);
                        }
                    }
                    else if (0 != C.coeffRef(i, j))
                        C.coeffRef(i, j) = 0;

                    if (ix <= kx)
                        ++it_i;
                    if (kx <= ix)
                        ++it_k;
                }

                for (vector<pair<int, int>>::iterator it = upd.begin(); it != upd.end(); ++it)
                    C.coeffRef(i, it->first) = it->second / local_gcd;
            }

            dividend = addInt(multiplyInt(Ckkr, b(i)), -multiplyInt(Cikr, b(k)));
            b(i) = divideInt(dividend, o);
        }
        if (!gcd_division)
        {
            o = Ckkr;
            divideRowByGCD(C, k);
        }
    }

    total_time = subtractTimes(start);
#ifndef TRULYSEAMLESS_SILENT
    cout << "running time: " << total_time << " ms" << endl;
#endif
}

/* Algorithm 2 IRREF by fraction-free Jordan Elimination */
void IRREF(SparseMatrixi& C, VectorXd& b, VectorXi indexC, bool gcd_division, double& total_time)
{
    auto start = std::chrono::high_resolution_clock::now();

    int dividend = 0;
    int o = 1;

    for (int k = C.outerSize() - 1; k > 0; k--)
    {
        // all zeroes row
        if (indexC(k) == -1)
            continue;

        // find pivot in row k
        int l = indexC(k);

        if (gcd_division)
            divideRowByGCD(C, k);

        // elimination(row i by j)
        int Ckl = C.coeff(k, l);

        for (int i = k - 1; i > -1; i--)
        {
            SparseMatrixi::InnerIterator it_i(C, i);
            int Cil = C.coeff(i, l);

            if (Cil == 0)
            {
                for (; it_i; ++it_i)
                {
                    if (it_i.index() < indexC(i) && it_i.value() != 0)
                        it_i.valueRef() = 0;

                    if (!gcd_division)
                    {
                        if (it_i.index() >= indexC(i) && it_i.value() != 0)
                        {
                            dividend = multiplyInt(Ckl, it_i.value());
                            it_i.valueRef() = divideInt(dividend, o);
                        }
                    }
                }
            }
            else
            {
                int local_gcd = -1;
                vector<pair<int, int>> upd;
                SparseMatrixi::InnerIterator it_k(C, k);
                while (it_i || it_k)
                {
                    while ((it_i) && it_i.value() == 0)
                        ++it_i;
                    while ((it_k) && it_k.value() == 0)
                        ++it_k;

                    if (!(it_i || it_k))
                        break;

                    int ix = (it_i) ? it_i.index() : C.cols();
                    int kx = (it_k) ? it_k.index() : C.cols();
                    int j = min(ix, kx);

                    if (j >= indexC(i))
                    {
                        int Cij = (ix == j) ? it_i.value() : 0;

                        int d = multiplyInt(Ckl, Cij);
                        if (j >= l)
                        {
                            int Ckj = (kx == j) ? it_k.value() : 0;
                            d = addInt(d, -multiplyInt(Cil, Ckj));
                        }
                        int v = divideInt(d, o);
                        upd.push_back(pair<int, int>(j, v));

                        if (v < 0)
                            v *= -1;
                        if (gcd_division && v > 0)
                        {
                            if (local_gcd < 0)
                                local_gcd = v;
                            else
                                local_gcd = gcd(local_gcd, v);
                        }
                    }

                    if (ix <= kx)
                        ++it_i;
                    if (kx <= ix)
                        ++it_k;
                }

                for (vector<pair<int, int>>::iterator it = upd.begin(); it != upd.end(); ++it)
                    C.coeffRef(i, it->first) = it->second / local_gcd;
            }

            dividend = addInt(multiplyInt(Ckl, b(i)), -multiplyInt(Cil, b(k)));
            b(i) = divideInt(dividend, o);
        }

        if (!gcd_division)
        {
            o = Ckl;
            divideRowByGCD(C, k);
        }
    }

    double t = subtractTimes(start);
    total_time += t;
#ifndef TRULYSEAMLESS_SILENT
    cout << "running time: " << t << " ms" << endl;
#endif
}


double getIntFactor(VectorXd& X_bar, double M, double uv_max)
{
    if (X_bar.size() <= 0)
        return 1.0;

    int delta = ceil(log2(uv_max));
    if (delta < 0)
        delta = 0;
    delta = 52 - delta - 1; // save on precision?
    if (delta <= 0)
        delta = 3;

    double x_max = fabs(X_bar(0));
    for (unsigned int i = 1; i < X_bar.size(); i++)
        if (x_max < fabs(X_bar(i)))
            x_max = fabs(X_bar(i));

    x_max = M / x_max;
    int log2_f = (int)floor(log2(x_max));

    return pow(2.0, log2_f);
}

/* Algorithm 3 Evaluation */
VectorXd
evaluate(SparseMatrixi C, VectorXd X_bar, VectorXd b, double M, double uv_max, VectorXi indexC, VectorXi indexR)
{
#ifndef TRULYSEAMLESS_SILENT
    auto start = std::chrono::high_resolution_clock::now();
#endif

    VectorXd X(C.innerSize());
    X.setZero();

    double f_max = getIntFactor(X_bar, M, uv_max);
    size_t max_S_size = 0;

    for (int k = C.innerSize() - 1; k > -1; k--)
    {
        if (indexR(k) == -1)
        {
            vector<int> D;
            D.resize(0);

            for (int i = 0; i < min((int)C.outerSize(), k + 1); i++)
                if (C.coeff(i, k) != 0 && indexC(i) != -1)
                    D.push_back(C.coeff(i, indexC(i)));

            // X(k) = makeDiv(fixedPrec(X_bar, k, f_max),D);
            X(k) = makeDiv(X_bar(k), D, f_max);
        }
        else
        {
            vector<pair<int, double>> S;
            S.resize(0);
            int j = indexR(k);
            int Cjk = C.coeff(j, k);

            for (SparseMatrixi::InnerIterator it(C, j); it; ++it)
            {
                if (it.index() > k && it.value() != 0 && X(it.index()) != 0)
                {
                    S.push_back(pair<int, double>(it.value(), X(it.index()) / Cjk));
                }
            }

            if (max_S_size < S.size())
                max_S_size = S.size();

            X(k) = b(j) - safeDotProd(S, M);
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    cout << "running time: " << subtractTimes(start) << " ms" << endl;
#endif
    return X / f_max;
}

/* Algorithm 4 safeDotProd */
double safeDotProd(vector<pair<int, double>> S, double M)
{
    stack<pair<int, double>> P;
    stack<pair<int, double>> N;

    for (size_t i = 0; i < S.size(); i++)
    {
        if (S[i].first * S[i].second > 0)
            P.push(pair<int, double>(std::abs(S[i].first), std::abs(S[i].second)));
        if (S[i].first * S[i].second < 0)
            N.push(pair<int, double>(std::abs(S[i].first), -std::abs(S[i].second)));
        if (S[i].first * S[i].second == 0)
        {
            cerr << "ERROR MESSAGE : zero element in set S!" << endl;
            cin.get();
        }
    }

    double r = 0;
    int counter = 1;
    while (!P.empty() || !N.empty())
    {
        counter++;
        if ((r < 0 && !P.empty()) || N.empty())
        {
            int c = P.top().first;
            double x = P.top().second;
            P.pop();

            int k = c < floor((M - r) / x) ? c : floor((M - r) / x);
            r += k * x;

            if (0 == counter % 1000 && 0 == k)
            {
#ifndef TRULYSEAMLESS_SILENT
                std::cout << "k = " << k << flush;
#endif
                throw std::runtime_error("Safe dot product numerical failure");
            }

            if (k < c)
                P.push(pair<int, double>(c - k, x));
        }
        else
        {
            int c = N.top().first;
            double x = N.top().second;
            N.pop();

            int k = c < floor((-M - r) / x) ? c : floor((-M - r) / x);
            r += k * x;
            if (0 == counter % 1000 && 0 == k)
            {
#ifndef TRULYSEAMLESS_SILENT
                std::cout << "k = " << k << flush;
#endif
                throw std::runtime_error("Safe dot product numerical failure");
            }

            if (k < c)
                N.push(pair<int, double>(c - k, x));
        }
    }

    return r;
}

/* Algorithm 5 makeDiv */
double makeDiv(double x, vector<int> D, double& f_max)
{
    int d = lcmArr(D);
    return fixedPrec(x / d, f_max) * d;
}

/* Algorithm 6 fixedPrec */
double fixedPrec(VectorXd& X_bar, int k, double& f_max)
{
    // return X_bar(k);
    double x = X_bar(k) * f_max;
    int x_int = x;
    return (double)x_int;
}

double fixedPrec(double x, double& f_max)
{
    x *= f_max;
    int x_int = x;
    return (double)x_int;
}

/* LCM */
int gcd(int a, int b)
{
    if (b)
        while ((a %= b) && (b %= a))
            ;
    return a + b;
}

int lcm(int a, int b)
{
    return multiplyInt(a, b) / gcd(a, b);
}

int lcmArr(vector<int> D)
{
    if (D.size() == 0)
        return 1;

    if (D.size() == 1)
        return D[0];

    int last_lcm = lcm(D[0], D[1]);
    for (size_t i = 2; i < D.size(); i++)
    {
        last_lcm = lcm(last_lcm, D[i]);
    }
    return last_lcm;
}

void statistic(SparseMatrixi C)
{
    int count = 0;
    int abs_max = 0;

    for (int i = 0; i < C.outerSize(); i++)
    {
        for (SparseMatrixi::InnerIterator it(C, i); it; ++it)
        {
            if (it.value() != 0)
            {
                count++;
                int v = it.value();
                if (v < 0)
                    v = v * (-1);
                if (v > abs_max)
                    abs_max = v;
            }
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    cout << "absolute largest value : " << abs_max << endl;
    cout << "entries:" << (double)100 * C.nonZeros() / (C.outerSize() * C.innerSize()) << "% : " << C.nonZeros()
         << endl;
    cout << "non-zero elements:" << (double)100 * count / (C.outerSize() * C.innerSize()) << "% : " << count << endl;
    cout << endl;
#endif
}

void statistic(SparseMatrixi C, VectorXi indexC)
{
    int count = 0;
    int abs_max = 0;

    for (int i = 0; i < C.outerSize(); i++)
    {
        for (SparseMatrixi::InnerIterator it(C, i); it; ++it)
        {
            if (it.index() >= indexC(i) && it.value() != 0)
            {
                count++;
                int v = it.value();
                if (v < 0)
                    v = v * (-1);
                if (v > abs_max)
                    abs_max = v;
            }
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    cout << "absolute largest value : " << abs_max << endl;
    cout << "entries:" << (double)100 * C.nonZeros() / (C.outerSize() * C.innerSize()) << "% : " << C.nonZeros()
         << endl;
    cout << "non-zero elements:" << (double)100 * count / (C.outerSize() * C.innerSize()) << "% : " << count << endl;
    cout << endl;
#endif
}

int divideInt(int dividend, int divisor, int s)
{
    (void)s;
    if (dividend % divisor != 0)
    {
#ifndef TRULYSEAMLESS_SILENT
        cerr << "ERROR MESSAGE div (" << s << ") : aliquant div!" << dividend << "/" << divisor << endl;
#endif
        throw std::runtime_error("Non-integer division result in divideInt");
    }

    return dividend / divisor;
}

int multiplyInt(int a, int b)
{
    if (a > 0)
    {
        if ((b > 0 && a > INT_MAX / b) || (b <= 0 && b < INT_MIN / a))
        {
#ifndef TRULYSEAMLESS_SILENT
            cerr << "ERROR MESSAGE : overflow Mult1!" << a << "x" << b << endl;
#endif
            throw std::runtime_error("Overflow in multiplyInt");
        }
    }
    else
    {
        if ((b > 0 && a < INT_MIN / b) || (b <= 0 && a != 0 && b < INT_MAX / a))
        {
#ifndef TRULYSEAMLESS_SILENT
            cerr << "ERROR MESSAGE : overflow Mult2!" << a << "x" << b << endl;
#endif
            throw std::runtime_error("Overflow in multiplyInt");
        }
    }
    return a * b;
}

int addInt(int a, int b)
{
    if ((b > 0 && a > INT_MAX - b) || (b < 0 && a < INT_MIN - b))
    {
#ifndef TRULYSEAMLESS_SILENT
        cerr << "ERROR MESSAGE : overflow Add!" << a << "+" << b << endl;
#endif
        throw std::runtime_error("Overflow in addInt");
    }
    return a + b;
}
