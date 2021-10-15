#ifndef TRAJECTORYCOSTFUNCTION_H
#define TRAJECTORYCOSTFUNCTION_H

#include <vector>
#include "ceres/ceres.h"

namespace SPL_OPTI
{
using ceres::CostFunction;
using ceres::Problem;
using ceres::SizedCostFunction;
using ceres::Solve;
using ceres::Solver;

class TimeCostFunction : public SizedCostFunction<1,1>
{
public:
  TimeCostFunction(double w)
  {
    m_w = w;
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {
    double beta = parameters[0][0];

    // Residule
    F[0] = m_w/beta;

    // Jacobian
    if (J != nullptr && J[0] != nullptr)
        J[0][0] = -m_w/(beta*beta);

    return true;
  }

private:
  double m_w;
};

class VelCostFunction : public SizedCostFunction<2,1,1,1>
{
public:
  VelCostFunction(double w, double w_violation, double max, double min)
  {
    m_w = w;
    m_w_violation = w_violation;
    m_max = max;
    m_min = min;
  }

  inline void assign_values(int i, double* F,
                            double** J,
                            const double &res,
                            const double &j0,
                            const double &j1,
                            const double &j2) const
  {
    // Residule
    F[i] = res;

    // Jacobian
    if (J != nullptr)
    {
      if (J[0] != nullptr)
        J[0][i] = j0;

      if (J[1] != nullptr)
        J[1][i] = j1;

      if (J[2] != nullptr)
        J[2][i] = j2;
    }
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {

    double q = parameters[1][0] - parameters[0][0];
    double beta = parameters[2][0];
    double c = m_w*beta;
    double c_violation = m_w_violation*beta;
    double vel = q*beta;

    assign_values(0,F,J,vel*m_w,-c,c,q*m_w);

    // Violation part
    if (vel > m_max)
    {
      assign_values(1,F,J,(vel - m_max)*m_w_violation,-c_violation,c_violation,q*m_w_violation);
    }
    else if (vel < m_min)
    {
      assign_values(1,F,J,(-vel + m_min)*m_w_violation,c_violation,-c_violation,-q*m_w_violation);
    }
    else
    {
      assign_values(1,F,J,0,0,0,0);
    }

    return true;
  }

private:
  double m_w, m_w_violation, m_max, m_min;
};

class AccCostFunction : public SizedCostFunction<2,1,1,1,1>
{
public:
  AccCostFunction(double w, double w_violation, double max, double min)
  {
    m_w = w;
    m_w_violation = w_violation;
    m_max = max;
    m_min = min;
  }

  inline void assign_values(int i, double* F,
                            double** J,
                            const double &res,
                            const double &j0,
                            const double &j1,
                            const double &j2,
                            const double &j3) const
  {
    // Residule
    F[i] = res;

    // Jacobian
    if (J != nullptr)
    {
      if (J[0] != nullptr)
        J[0][i] = j0;

      if (J[1] != nullptr)
        J[1][i] = j1;

      if (J[2] != nullptr)
        J[2][i] = j2;

      if (J[3] != nullptr)
        J[3][i] = j3;
    }
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {

    double q = parameters[2][0] - 2*parameters[1][0] + parameters[0][0];
    double beta = parameters[3][0];
    double c = m_w*beta*beta;
    double c_violation = m_w_violation*beta*beta;
    double acc = q*beta*beta;

    assign_values(0,F,J,acc*m_w,c,-2*c,c,2*m_w*beta*q);

    // Violation part
    if (acc > m_max)
    {
      assign_values(1,F,J,(acc - m_max)*m_w_violation,
                    c_violation,-2*c_violation,c_violation,2*m_w_violation*beta*q);
    }
    else if (acc < m_min)
    {
      assign_values(1,F,J,(-acc + m_min)*m_w_violation,
                    -c_violation,2*c_violation,-c_violation,-2*m_w_violation*beta*q);
    }
    else
    {
      assign_values(1,F,J,0,0,0,0,0);
    }
    return true;
  }

private:
  double m_w, m_w_violation, m_max, m_min;
};

class JerkCostFunction : public SizedCostFunction<2,1,1,1,1,1>
{
public:
  JerkCostFunction(double w, double w_violation, double max, double min)
  {
    m_w = w;
    m_w_violation = w_violation;
    m_max = max;
    m_min = min;
  }

  inline void assign_values(int i, double* F,
                            double** J,
                            const double &res,
                            const double &j0,
                            const double &j1,
                            const double &j2,
                            const double &j3,
                            const double &j4) const
  {
    // Residule
    F[i] = res;

    // Jacobian
    if (J != nullptr)
    {
      if (J[0] != nullptr)
        J[0][i] = j0;

      if (J[1] != nullptr)
        J[1][i] = j1;

      if (J[2] != nullptr)
        J[2][i] = j2;

      if (J[3] != nullptr)
        J[3][i] = j3;

      if (J[4] != nullptr)
        J[4][i] = j4;
    }
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {

    double q = parameters[3][0] - 3*parameters[2][0] + 3*parameters[1][0] - parameters[0][0];
    double beta = parameters[4][0];
    double c = m_w*beta*beta*beta;
    double c_violation = m_w_violation*beta*beta*beta;
    double jerk = q*beta*beta*beta;

    assign_values(0,F,J,jerk*m_w,-c,3*c,-3*c,c,3*m_w*beta*beta*q);

    // Violation part
    if (jerk > m_max)
    {
      assign_values(1,F,J,(jerk - m_max)*m_w_violation,
                    -c_violation,3*c_violation,-3*c_violation,c_violation,3*m_w_violation*beta*beta*q);
    }
    else if (jerk < m_min)
    {
      assign_values(1,F,J,(-jerk + m_min)*m_w_violation,
                    c_violation,-3*c_violation,3*c_violation,-c_violation,-3*m_w_violation*beta*beta*q);
    }
    else
    {
      assign_values(1,F,J,0,0,0,0,0,0);
    }
    return true;
  }

private:
  double m_w, m_w_violation, m_max, m_min;
};

class StartCostFunction : public SizedCostFunction<3,1,1,1,1>
{
public:
  StartCostFunction(double w)
  {
    m_w = w;
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {
    Eigen::Vector3d cp(parameters[0][0],parameters[1][0],parameters[2][0]);
    double beta = parameters[3][0];

    Eigen::Matrix3d M;
    M<<1.0/6.0, 2.0/3.0, 1.0/6.0,
        -0.5*beta, 0, 0.5*beta,
        beta*beta, -2*beta*beta, beta*beta;

    Eigen::Vector3d inip = M*cp;

    // Residule
    for(int i=0; i<3; i++)
      F[i] = m_w*(inip(i)-m_s_ini(i));

    // Jacobian
    if (J != nullptr)
    {
      if (J[0] != nullptr)
      {
        J[0][0] = m_w/6;
        J[0][1] = -beta*m_w/2;
        J[0][2] = beta*beta*m_w;
      }

      if (J[1] != nullptr)
      {
        J[1][0] = 2*m_w/3;
        J[1][1] = 0;
        J[1][2] = -2*beta*beta*m_w;
      }

      if (J[2] != nullptr)
      {
        J[2][0] = m_w/6;
        J[2][1] = beta*m_w/2;
        J[2][2] = beta*beta*m_w;
      }

      if (J[3] != nullptr)
      {
        J[3][0] = 0;
        J[3][1] = -m_w*(cp(0) - cp(2))/2;
        J[3][2] = 2*beta*m_w*(cp(0) - 2*cp(1) + cp(2));
      }
    }
    return true;
  }

public:
  Eigen::Vector3d m_s_ini;

private:
  double m_w;
};

class FinishCostFunction : public SizedCostFunction<3,1,1,1,1>
{
public:
  FinishCostFunction(double w)
  {
    m_w = w;
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* F,
                        double** J) const
  {
    Eigen::Vector3d cp(parameters[0][0],parameters[1][0],parameters[2][0]);
    double beta = parameters[3][0];

    Eigen::Matrix3d M;
    M<<1.0/6.0, 2.0/3.0, 1.0/6.0,
        -0.5*beta, 0, 0.5*beta,
        beta*beta, -2*beta*beta, beta*beta;

    Eigen::Vector3d terp = M*cp;

    // Residule
    for(int i=0; i<3; i++)
      F[i] = m_w*(terp(i)-m_s_ter(i));

    // Jacobian
    if (J != nullptr)
    {
      if (J[0] != nullptr)
      {
        J[0][0] = m_w/6;
        J[0][1] = -beta*m_w/2;
        J[0][2] = beta*beta*m_w;
      }

      if (J[1] != nullptr)
      {
        J[1][0] = 2*m_w/3;
        J[1][1] = 0;
        J[1][2] = -2*beta*beta*m_w;
      }

      if (J[2] != nullptr)
      {
        J[2][0] = m_w/6;
        J[2][1] = beta*m_w/2;
        J[2][2] = beta*beta*m_w;
      }

      if (J[3] != nullptr)
      {
        J[3][0] = 0;
        J[3][1] = -m_w*(cp(0) - cp(2))/2;
        J[3][2] = 2*beta*m_w*(cp(0) - 2*cp(1) + cp(2));
      }
    }
    return true;
  }

public:
  Eigen::Vector3d m_s_ter;

private:
  double m_w;
};
}

#endif // TRAJECTORYCOSTFUNCTION_H
