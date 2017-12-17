#include "MVMO.h"
#include <iostream>
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

int main()
{
    VectorXd lb(2), ub(2);
    lb << -1, -1;
    ub << 1, 1;
    MVMO optimizer(quad, lb, ub);
    optimizer.optimize();
    cout << optimizer.best_x().transpose() << endl;
    cout << optimizer.best_y() << endl;
    return EXIT_SUCCESS;
}
