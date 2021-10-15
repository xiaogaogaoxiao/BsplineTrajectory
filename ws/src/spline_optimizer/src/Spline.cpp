#include "spline_optimizer/Spline.h"
namespace SPL_OPTI
{
Spline::Spline()
{

}

void Spline::init(int p, int no_of_ctrl_pnts, double beta, int D, bool optimize_beta)
{
  m_p = p;
  m_n = no_of_ctrl_pnts - 1;
  m_beta = beta;
  m_m = m_p + m_n + 1;
  m_u.resize(m_m+1);
  for (int i=0; i<=m_m; i++)
    m_u[i] = i;
  m_D = D;
  m_optimize_beta = optimize_beta;

  if (m_p == 3)
  {
    calc_Q_v();
    calc_Q_a();
    calc_Q_j();
  }

  m_v_max.resize(m_D);
  m_v_max<<1,1,1;

  m_v_min.resize(m_D);
  m_v_min<<-1,-1,-1;

  m_a_max.resize(m_D);
  m_a_max<<1,1,1;

  m_a_min.resize(m_D);
  m_a_min<<-1,-1,-1;

  m_j_max.resize(m_D);
  m_j_max<<1,1,1;

  m_j_min.resize(m_D);
  m_j_min<<-1,-1,-1;
}

void Spline::set_control_points(const Eigen::MatrixXd &ctrl_pnts)
{
  if (ctrl_pnts.rows() != m_n + 1)
    std::wcerr<<"The size of the control point is wrong!"<<std::endl;

  if (ctrl_pnts.cols() != m_D)
    std::wcerr<<"The dimension of the control point is wrong!"<<std::endl;

  m_ctrl_points = ctrl_pnts;
}

Eigen::VectorXd Spline::single_deboor(double u_probe)
{
  // Bound the u_probe
  u_probe = std::min(std::max(m_u(m_p), u_probe), m_u(m_m - m_p));

  // find which knot segment u_probe belongs to
  int k = m_p;

  while (true)
  {
    if (k+1 > m_u.size())
      std::wcerr<<"Out side the knot range!"<<std::endl;

    if (m_u(k+1) >= u_probe)
      break;

    k++;
  }

  // t(t_ctt) is maped to knot segment u_k ~ u_{k+1}
  // there are at most p+1 basis functions N_k-p,p(u), N_k-p+1,p(u), ..., N_k,p(u) non-zero on knot span [uk, uk+1)
  // the effective control points are P_k-p ~ P_k
  // since Matlab index start from 1 instead of 0
  // The effective control points are
  Eigen::MatrixXd d = m_ctrl_points.block(k - m_p, 0, m_p+1, m_D);

  for (int r = 1; r <= m_p; r++)
  {
    for (int i= m_p; i >= r; i--)
    {
      double alpha = (u_probe - m_u(i + k - m_p)) / (m_u(i + 1 + k - r) - m_u(i + k - m_p));
      d.row(i) = (1-alpha) * d.row(i-1) + alpha * d.row(i);
    }
  }

  return d.row(m_p);
}

Eigen::MatrixXd Spline::get_trajectory(std::vector<double> t)
{
  Eigen::MatrixXd trajectory;
  trajectory.resize(t.size(),m_D);

  for (size_t i = 0; i < t.size(); i++)
  {
    double u_probe = t[i] * m_beta + m_u(m_p);
    trajectory.row(i) = single_deboor(u_probe);
  }
  return trajectory;
}

Eigen::Vector2d Spline::get_available_s_range()
{
  Eigen::Vector2d range;
  range << m_u(m_p), m_u(m_m-m_p);
  return range;
}

Eigen::Vector2d Spline::get_available_t_range()
{
  Eigen::Vector2d range;
  range << 0, m_u(m_m-m_p)-m_u(m_p);
  return range/m_beta;
}

double Spline::N_func(int j, int k, double s)
{
  if (j==0 && k == 0)
    return 1.0;
  else if (j == 0 && k != 0)
    return (1.0-s)/k*N_func(0,k-1,s);
  else if (j == k && k != 0)
    return s/k*N_func(j-1,k-1,s);
  else
    return (k-j+s)/k*N_func(j-1,k-1,s) + (1.0+j-s)/k*N_func(j,k-1,s);
}

double Spline::B_func(int k, double s)
{
  double s_lb = floor(s);
  if (s_lb >=0 && s_lb <= k)
    return N_func(static_cast<int>(k-s_lb), k, s-s_lb);
  else
    return 0;
}

void Spline::construct_evaluation_matrix(std::vector<double> s)
{
  m_Eva.resize(s.size(), m_n + 1);

  for(int i=0; i<s.size(); i++)
  {
    for (int j=0; j<m_n+1; j++)
    {
      m_Eva(i,j) = B_func(m_p, s[i] - m_u(j));
    }
  }
}

void Spline::construct_evaluation_matrix(const Eigen::VectorXd &s)
{
  std::vector<double> s_tilt;
  for(int i=0; i<s.size(); i++)
    s_tilt.push_back(s[i]);

  construct_evaluation_matrix(s_tilt);
}

void Spline::construct_evaluation_matrix(int n)
{
  std::vector<double> s;
  Eigen::VectorXd s_tilt;
  s_tilt.resize(n);
  s_tilt.setLinSpaced(n,m_u(m_p),m_u(m_m-m_p));

  for(int i=0;i<n;i++)
    s.push_back(s_tilt[i]);

  construct_evaluation_matrix(s);
}

Spline Spline::get_derivative()
{
  Spline dS = *this;
  dS.m_p -= 1;
  dS.m_n -= 1;
  dS.m_m = dS.m_p + dS.m_n + 1;
  // remove the head and tail (m_u originally have m_m+1 entries)
  dS.m_u = m_u.block(1,0,m_m-1,1);
  dS.m_ctrl_points.resize(dS.m_n + 1, m_D);
  for (int i=0; i<m_n; i++)
  {
    dS.m_ctrl_points.row(i) = m_beta*(m_ctrl_points.row(i+1) - m_ctrl_points.row(i));
  }
  return dS;
}

void Spline::set_ini_ter_matrix()
{
  m_A_ini<<1.0/6.0, 2.0/3.0, 1.0/6.0,
           -0.5*m_beta, 0, 0.5*m_beta,
           m_beta*m_beta, -2*m_beta*m_beta, m_beta*m_beta;

  m_A_ter<<1.0/6.0, 2.0/3.0, 1.0/6.0,
           -0.5*m_beta, 0, 0.5*m_beta,
           m_beta*m_beta, -2*m_beta*m_beta, m_beta*m_beta;

}

void Spline::init_with_approximation(const Eigen::MatrixXd &s_ini, const Eigen::MatrixXd &s_ter,
                             const Eigen::MatrixXd &D, const Eigen::VectorXd &S)
{
  if (m_p != 3)
  {
    std::wcerr<<"This function requires the spline has order 3!"<<std::endl;
    return;
  }

  set_ini_ter_matrix();

  // first determine the first and last three control points from the
  // initial and terminal conditions
  m_ctrl_points.setZero(m_n+1,m_D);
  Eigen::VectorXi fix;
  fix.setZero(m_n+1);

  m_ctrl_points.block(0,0,3,m_D) = m_A_ini.colPivHouseholderQr().solve(s_ini);
  fix.block(0,0,3,1) << 1,1,1;

  m_s_ini = s_ini;
  m_s_ter = s_ter;

  //similarly for s_ter
  m_ctrl_points.block(m_n-2,0,3,m_D) = m_A_ter.colPivHouseholderQr().solve(s_ter);
  fix.block(m_n-2,0,3,1) << 1,1,1;

  // execute the actual approximation
  approximation(D,S,fix);
}

void Spline::approximation(const Eigen::MatrixXd &D, const Eigen::VectorXd &S, const Eigen::MatrixXi &fix)
{
  // Construct the mapping matrix from fix
  int k=0;
  Eigen::MatrixXd M;
  M.setZero(m_n+1, m_n+1);

  for (int i=0; i<m_n+1; i++)
  {
    if (fix(i) == 1)
    {
      M(k,i) = 1;
      k++;
    }
  }
  int k_fix = k;
  // Get C_f for later use (C_f is the fixed the control points)
  Eigen::MatrixXd C_f = (M*m_ctrl_points).block(0,0,k_fix,m_D);

  // The second scan
  for (int i=0; i<m_n+1; i++)
  {
    if (fix(i) == 0)
    {
      M(k,i) = 1;
      k++;
    }
  }
  M.transposeInPlace();

  //Construct the rest of the matrices
  construct_evaluation_matrix(S);
  Eigen::MatrixXd H = m_Eva.transpose()*m_Eva + 0.5*m_beta*m_Q_v + 0.5*m_beta*m_beta*m_beta*m_Q_a;
  Eigen::MatrixXd F = -2*D.transpose()*m_Eva;
  H = M.transpose()*H*M;
  F = F*M;
  Eigen::MatrixXd H_fv = H.block(0, k_fix, k_fix, m_n-k_fix+1);
  Eigen::MatrixXd H_vv = H.block(k_fix,k_fix, m_n-k_fix+1, m_n-k_fix+1);
  Eigen::MatrixXd F_v  = F.block(0,k_fix,m_D,m_n-k_fix+1);
  Eigen::MatrixXd C_v  = H_vv.colPivHouseholderQr().solve(-0.5*(2*C_f.transpose()*H_fv+F_v).transpose());
  Eigen::MatrixXd C(C_f.rows()+C_v.rows(), C_f.cols());
  C<<C_f,C_v;
  m_ctrl_points = M*C;

//  std::cout<<m_ctrl_points<<std::endl;
}

void Spline::get_params_from_array(const double *x)
{
  int k = 0;
  for (int d=0; d<m_D; d++)
  {
    for (int i=0; i<m_n+1; i++)
    {
      m_ctrl_points(i,d) = x[k++];
    }
  }

  if (m_optimize_beta)
  {
    m_beta = x[k];
  }
}

void Spline::send_params_to_array(double *x)
{
  int k = 0;
  for (int d=0; d<m_D; d++)
  {
    for (int i=0; i<m_n+1; i++)
    {
      x[k++] = m_ctrl_points(i,d);
    }
  }

  if (m_optimize_beta)
  {
    x[k] = m_beta;
  }
}

int Spline::total_parameter_size()
{
  // If m_optimize_beta flag is set
  // the optimization parameter now contains m_beta
  // therefore add one here.
  if (m_optimize_beta)
    return m_D*(m_n+1) + 1;
  else
    return m_D*(m_n+1);
}

// i is the id of the ctrl point, d is its dimension (both start
// from 1)
int Spline::get_rowid(int i, int d)
{
  return i*m_D + d;
}
// i is the id of the ctrl point, d is its dimension (both start
// from 1)
int Spline::get_colid(int i, int d)
{
  return d*(m_n+1) + i;
}
//------------------------------------------------------------//
//                       COST FUNCTIONS                       //
//------------------------------------------------------------//
void Spline::add_vel_cost(ceres::Problem *problem)
{
  for (int i=0;i<m_n;i++)
  {
    for(int d=0; d<m_D; d++)
    {
      problem->AddResidualBlock(new SPL_OPTI::VelCostFunction(0.5,5,1,-1),nullptr,
                               &m_ctrl_points(i,d), &m_ctrl_points(i+1,d), &m_beta);
    }
  }
}

void Spline::add_acc_cost(ceres::Problem* problem)
{
  for (int i=0;i<m_n-1;i++)
  {
    for(int d=0; d<m_D; d++)
    {
      problem->AddResidualBlock(new SPL_OPTI::AccCostFunction(0.5,5,2,-2),nullptr,
                               &m_ctrl_points(i,d), &m_ctrl_points(i+1,d), &m_ctrl_points(i+2,d), &m_beta);
    }
  }
}

void Spline::add_jerk_cost(ceres::Problem* problem)
{
  for (int i=0;i<m_n-2;i++)
  {
    for(int d=0; d<m_D; d++)
    {
      problem->AddResidualBlock(new SPL_OPTI::JerkCostFunction(0.5,5,2,-2),nullptr,
                               &m_ctrl_points(i,d), &m_ctrl_points(i+1,d), &m_ctrl_points(i+2,d),
                               &m_ctrl_points(i+3,d), &m_beta);
    }
  }
}

void Spline::add_start_cost(ceres::Problem* problem)
{
  for (int d=0; d<m_D; d++)
  {
    SPL_OPTI::StartCostFunction *cf_start= new SPL_OPTI::StartCostFunction(30);
    cf_start->m_s_ini = m_s_ini.col(d);
    problem->AddResidualBlock(cf_start,nullptr,
                             &m_ctrl_points(0,d), &m_ctrl_points(1,d), &m_ctrl_points(2,d),&m_beta);
  }
}

void Spline::add_finish_cost(ceres::Problem* problem)
{
  for (int d=0; d<m_D; d++)
  {
    SPL_OPTI::FinishCostFunction *cf_finish= new SPL_OPTI::FinishCostFunction(30);
    cf_finish->m_s_ter = m_s_ter.col(d);
    problem->AddResidualBlock(cf_finish,nullptr,
                             &m_ctrl_points(m_n-2,d), &m_ctrl_points(m_n-1,d), &m_ctrl_points(m_n,d),&m_beta);
  }
}

void Spline::add_time_cost(ceres::Problem *problem)
{
  problem->AddResidualBlock(new SPL_OPTI::TimeCostFunction(20.0),nullptr,&m_beta);
}

}
