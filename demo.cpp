#include "MVMO.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace Eigen;

double quad(const VectorXd& xs)
{
    double y = 0.0;
    for (long i = 0; i < xs.size(); ++i)
    {
        y += pow(xs(i), 2);
    }
    return y;
}
double rosenbrock(const VectorXd& xs)
{
    double y = 0.0;
    assert(xs.size() == 2);
    const double a  = 1;
    const double b  = 100;
    const double x1 = xs(0);
    const double x2 = xs(1);
    return pow(a - x1, 2) + b * pow(x2 - x1 * x1, 2);
}

double ackley(const VectorXd& xs_)
{
    VectorXd xs = xs_;
    for(long i = 0; i < xs.size(); ++i)
        xs(i) -= i;
    const double a = 20;
    const double b = 0.2;
    const double c = 2 * M_PI;
    double sum1 = 0, sum2 = 0;
    for(long i = 0; i < xs.size(); ++i)
    {
        const double xi = xs(i);
        sum1 += pow(xi, 2);
        sum2 += cos(c * xi);
    }
    const double term1 = -1 * a * exp(-b * sqrt(sum1 / xs.size()));
    const double term2 = -1 * exp(sum2 / xs.size());
    const double y     = term1 + term2 + a + exp(1);
    return y;
}

// double cec14(const VectorXd& xs)
// {
//     ofstream param("./circuit/param");
//     for(long i = 1; i <= xs.size(); ++i)
//     {
//         param << setprecision(18) << ".param x" << i << " = " << xs[i-1] << endl;
//     }
//     param.close();
//     system("cd circuit && perl run.pl");

//     ifstream result("./circuit/result.po");
//     double fom;
//     result >> fom;
//     result.close();
//     return fom;
// }

int main(int args_num, char** args)
{
    size_t dim      = 6;
    size_t num_eval = dim * 50;
    cout << setprecision(18);
    if(args_num > 1)
        num_eval = atoi(args[1]);
    const VectorXd lb = VectorXd::Constant(dim, 1, -32);
    const VectorXd ub = VectorXd::Constant(dim, 1, 32);
    MVMO optimizer(ackley, lb, ub);
    // these settings are optional
    optimizer.set_archive_size(25);
    optimizer.set_max_eval(num_eval);
    optimizer.set_fs_init(0.5);
    optimizer.set_fs_final(20);
    optimizer.optimize();
    cout << "BestX: " << optimizer.best_x().transpose() << endl;
    cout << "BestY: " << optimizer.best_y() << endl;
    return EXIT_SUCCESS;
}
