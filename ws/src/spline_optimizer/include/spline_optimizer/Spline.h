#ifndef Spline_H
#define Spline_H
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <ceres/ceres.h>
#include <spline_optimizer/TrajectoryCostFunction.h>

namespace SPL_OPTI
{
class Spline
{

public:
  Spline();
  void init(int p, int no_of_ctrl_pnts, double beta, int D, bool optimize_beta = false);
  void set_control_points(const Eigen::MatrixXd &ctrl_pnts);
  Eigen::VectorXd single_deboor(double u_probe);
  Eigen::MatrixXd get_trajectory(std::vector<double> t);
  Eigen::Vector2d get_available_s_range();
  Eigen::Vector2d get_available_t_range();
  void construct_evaluation_matrix(std::vector<double> s);
  void construct_evaluation_matrix(const Eigen::VectorXd &s);
  void construct_evaluation_matrix(int n);
  Spline get_derivative();
  void set_ini_ter_matrix();
  void init_with_approximation(const Eigen::MatrixXd &s_ini, const Eigen::MatrixXd &s_ter,
                               const Eigen::MatrixXd &D, const Eigen::VectorXd &S);
  void approximation(const Eigen::MatrixXd &D, const Eigen::VectorXd &S, const Eigen::MatrixXi &fix);
  void get_params_from_array(const double *x);
  void send_params_to_array(double*x);
  int total_parameter_size();
  int get_param_idx(int i, int d)
  {
    return get_colid(i,d);
  }

  void add_vel_cost(ceres::Problem* problem);
  void add_acc_cost(ceres::Problem* problem);
  void add_jerk_cost(ceres::Problem* problem);
  void add_start_cost(ceres::Problem* problem);
  void add_finish_cost(ceres::Problem* problem);
  void add_time_cost(ceres::Problem* problem);

private:
  double N_func(int j, int k, double s);
  double B_func(int k, double s);
  int get_rowid(int i, int d);
  int get_colid(int i, int d);

  void calc_Q_v()
  {
    int nm = m_n + 1;
    m_Q_v.resize(nm, nm);
    m_Q_v.setZero();

    if (m_p != 3)
    {
      std::wcerr<<"This function requires the spline has order 3!"<<std::endl;
      return;
    }

    for (int i = 0; i < nm; ++i)
    {
      if (i == 0)
        m_Q_v.block(i,i,1,4) << 6, 7,-12,-1;
      else if (i == 1)
        m_Q_v.block(i,i-1,1,5)<<7, 40, -22, -24, -1;
      else if (i == 2)
        m_Q_v.block(i,i-2,1,6)<<-12, -22, 74, -15, -24, -1;
      else if (i == nm-1)
        m_Q_v.block(i,i-3,1,4)<<-1, -12, 7, 6;
      else if (i == nm-2)
        m_Q_v.block(i,i-3,1,5)<<-1, -24, -22, 40, 7;
      else if (i == nm-3)
        m_Q_v.block(i,i-3,1,6)<<-1, -24, -15, 74, -22, -12;
      else
        m_Q_v.block(i,i-3,1,7)<<-1, -24, -15, 80, -15, -24, -1;
    }
    m_Q_v = m_Q_v/120.0;
  }

  void calc_Q_a()
  {
    int nm = m_n + 1;
    m_Q_a.resize(nm, nm);
    m_Q_a.setZero();

    if (m_p != 3)
    {
      std::wcerr<<"This function requires the spline has order 3!"<<std::endl;
      return;
    }

    for (int i = 0; i < nm; ++i)
    {
      if (i == 0)
        m_Q_a.block(i,i,1,4) << 2, -3, 0, 1;
      else if (i == 1)
        m_Q_a.block(i,i-1,1,5)<<-3, 8, -6, 0, 1;
      else if (i == 2)
        m_Q_a.block(i,i-2,1,6)<<0, -6, 14, -9, 0, 1;
      else if (i == nm-1)
        m_Q_a.block(i,i-3,1,4)<<1, 0, -3, 2;
      else if (i == nm-2)
        m_Q_a.block(i,i-3,1,5)<<1, 0, -6, 8, -3;
      else if (i == nm-3)
        m_Q_a.block(i,i-3,1,6)<<1, 0, -9, 14, -6, 0;
      else
        m_Q_a.block(i,i-3,1,7)<<1, 0, -9, 16, -9, 0, 1;
    }
    m_Q_a = m_Q_a/6.0;
  }

  void calc_Q_j()
  {
    int nm = m_n + 1;
    m_Q_j.resize(nm, nm);
    m_Q_j.setZero();

    if (m_p != 3)
    {
      std::wcerr<<"This function requires the spline has order 3!"<<std::endl;
      return;
    }

    for (int i = 0; i < nm; ++i)
    {
      if (i == 0)
        m_Q_j.block(i,i,1,4) << 1, -3, 3, -1;
      else if (i == 1)
        m_Q_j.block(i,i-1,1,5)<<-3, 10, -12, 6, -1;
      else if (i == 2)
        m_Q_j.block(i,i-2,1,6)<<3, -12, 19, -15, 6, -1;
      else if (i == nm-1)
        m_Q_j.block(i,i-3,1,4)<<-1, 3, -3, 1;
      else if (i == nm-2)
        m_Q_j.block(i,i-3,1,5)<<-1, 6, -12, 10, -3;
      else if (i == nm-3)
        m_Q_j.block(i,i-3,1,6)<<-1, 6, -15, 19, -12, 3;
      else
        m_Q_j.block(i,i-3,1,7)<<-1, 6, -15, 20, -15, 6, -1;
    }
  }

public:
  int m_p; // Degree
  int m_n; // n+1 is the number of control points
  int m_m; // m+1 is the number of knots, there is m = n + p +1;
  Eigen::VectorXd m_u; // The uniform knot vector [u0, u1, ..., um], u_{i+1} - u_i = 1, the valid domain is [up u_{m-p}]
  double m_beta; //  Time scale t(real time) * beta = u
  Eigen::MatrixXd m_ctrl_points ;// In n+1 x D format
  int m_D; // Dimension of the control points
  Eigen::MatrixXd m_Eva; // Evaluation matrix
  bool m_optimize_beta; // A flag to set whether to optimize beta, default true

  // The following only works for the case p = 3
  Eigen::MatrixXd m_Q_v; // beta*ctrl_points'*Q_v*ctrl_points is the /int v^2 dt
  Eigen::MatrixXd m_Q_a; // beta^3*ctrl_points'*Q_v*ctrl_points is the /int a^2 dt
  Eigen::MatrixXd m_Q_j; // beta^5*ctrl_points'*Q_v*ctrl_points is the /int j^2 dt
  Eigen::Matrix3d m_A_ini; // initial condition matrix A_ini*[cp1,cp2,cp3]' = [p0 v0 a0]
  Eigen::Matrix3d m_A_ter; // terminal condition matrix
  Eigen::MatrixXd m_s_ini; // initial state in the following form [px py pz; vx vy vz; ax ay az;]
  Eigen::MatrixXd m_s_ter; // terminal state in the following form [px py pz; vx vy vz; ax ay az;]

  // Velocity, acceleration and jerk limit
  Eigen::VectorXd m_v_max, m_v_min, m_a_max, m_a_min, m_j_max, m_j_min;
};
}
#endif // Spline_H
