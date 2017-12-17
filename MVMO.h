#pragma once
#include <Eigen/Dense>
#include <vector>
class MVMO
{
public:
    // TODO: Support constrained optimization
    typedef std::function<double(const Eigen::VectorXd&)> MVMO_Obj;
    MVMO(MVMO_Obj f, const Eigen::VectorXd& lb, const Eigen::VectorXd& ub);
    void optimize();
    void optimize_one_step(); 
    Eigen::VectorXcd best_x() const;
    double best_y() const;

    // set algorithm parameters
    void set_max_eval(size_t v) { _max_eval = v; }
    void set_num_init(size_t v) { _num_init = v; }

protected:
    const size_t _dim;
    const Eigen::VectorXd& _lb;
    const Eigen::VectorXd& _ub;
    size_t _eval_counter;

    size_t _max_eval;
    size_t _num_init;

    Eigen::VectorXd& _scale(const Eigen::VectorXd&) const;       // scale from [lb, ub] to [0, 1]
    Eigen::VectorXd& _scale_back(const Eigen::VectorXd&) const;  // scale from [0, 1] to [lb, ub]
    Eigen::VectorXd _convert(std::vector<double>&) const;
    std::vector<double> _convert(Eigen::VectorXd&) const;

    double _run_func(const Eigen::VectorXd& x);            // wrapper of the original objective function

private:
    MVMO_Obj _f;
};
