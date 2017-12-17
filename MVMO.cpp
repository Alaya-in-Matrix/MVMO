#include "MVMO.h"
#include "MVMO_util.h"
using namespace std;
using namespace Eigen;

MVMO::MVMO(MVMO::MVMO_Obj f, const VectorXd& lb, const VectorXd& ub)
    : _dim(lb.size()), _lb(lb), _ub(ub)
{}
void MVMO::optimize(){}
VectorXcd MVMO::best_x() const {}
double MVMO::best_y() const {}
