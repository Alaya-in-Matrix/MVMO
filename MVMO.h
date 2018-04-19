#pragma once
#include <Eigen/Dense>
#include <vector>
#include <random>
class MVMO
{
public:
    // TODO: Support constrained optimization
    typedef std::function<double(const Eigen::VectorXd&)> MVMO_Obj;
    MVMO(MVMO_Obj f, const Eigen::VectorXd& lb, const Eigen::VectorXd& ub);
    void optimize(const Eigen::MatrixXd& guess);
    void optimize();
    void optimize_one_step(); 
    Eigen::VectorXd best_x() const;
    double best_y() const;
    void info_print() const;

    // set algorithm parameters
    void set_max_eval(size_t v) { _max_eval = v; }
    void set_num_init(size_t v) { _num_init = v; }
    void set_fs_init(double f) { _fs_init = f; }
    void set_fs_final(double f) { _fs_final = f; }
    void set_m_init(size_t m) { _m_init = m; }
    void set_m_final(size_t m) { _m_final = m; }
    void set_archive_size(size_t as) { _archive_size = as; }
    void set_delta_d0(double v) { _delta_d0 = v; }
protected:

    // algorithm settings andparameters
    const size_t _dim;
    const Eigen::VectorXd& _lb;
    const Eigen::VectorXd& _ub;
    size_t _max_eval;
    size_t _num_init;
    size_t _archive_size;
    double _fs_init; // shape factor
    double _fs_final;
    size_t _m_init;
    size_t _m_final;
    double _delta_d0;

    // inner state variables
    size_t _eval_counter = 0;
    double _fshape;
    Eigen::MatrixXd _dbx; // dim * _eval_counter
    Eigen::VectorXd _dby; // _eval_counter * 1
    Eigen::MatrixXd _archive_x;
    Eigen::VectorXd _archive_mean;
    Eigen::VectorXd _archive_s;
    Eigen::VectorXd _archive_d;
    Eigen::VectorXd _archive_s1;
    Eigen::VectorXd _archive_s2;

#ifdef DEBUG_RAND_SEED
    std::mt19937_64 _engine = std::mt19937_64(DEBUG_RAND_SEED);
#else
    std::mt19937_64 _engine = std::mt19937_64(std::random_device{}());
#endif

    Eigen::VectorXd _best_x; //
    double _best_y = std::numeric_limits<double>::infinity(); //
    
    Eigen::VectorXd _scale(const Eigen::VectorXd&) const;       // scale from [lb, ub] to [0, 1]
    Eigen::VectorXd _scale_back(const Eigen::VectorXd&) const;  // scale from [0, 1] to [lb, ub]
    Eigen::VectorXd     _convert(std::vector<double>&) const;
    std::vector<double> _convert(Eigen::VectorXd&) const;

    double _run_func(const Eigen::VectorXd& x);            // wrapper of the original objective function
    Eigen:: VectorXd _run_func_batch(const Eigen::MatrixXd& xs);            // wrapper of the original objective function
    void _default_setting();
    void _init_archive();
    void _initialize();
    void _initialize(const Eigen::MatrixXd& guess);
    void   _fs();
    double _m();
    std::vector<size_t> _seq_idx(size_t) const;
    Eigen::MatrixXd _slice_matrix(const Eigen::MatrixXd&, const std::vector<size_t>&) const;
    void _archive();
    std::vector<size_t> _pick_from_seq(size_t, size_t);
    double _rand01();
    double _hfunc(double xbar, double s1, double s2, double x) const;
    Eigen::Vector2d _mean_var_noeq(const Eigen::RowVectorXd& xs) const;

private:
    MVMO_Obj _f;
};
