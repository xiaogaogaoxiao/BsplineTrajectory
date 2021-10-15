#include "ros/ros.h"
#include "geometry_msgs/Twist.h"
#include <fstream>
#include <deque>
#include <utility>
#include <nav_msgs/Odometry.h>
#include "tf/tf.h"
#include <sensor_msgs/Imu.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Bool.h>
#include "spline_optimizer/Spline.h"
#include "ceres/ceres.h"
#include "spline_optimizer/TrajectoryCostFunction.h"
#include "chrono"

int main(int argc, char **argv)
{
  ros::init(argc, argv, "spline_optimizer");
  ros::NodeHandle nh;

  SPL_OPTI::Spline S;
  int dimension = 3;
  S.init(3,40,1.0,dimension,true);

  Eigen::MatrixXd s_ini,s_ter,D;
  Eigen::VectorXd DS;
  s_ini.resize(3,dimension);
  s_ter.resize(3,dimension);
  s_ini<<0,0,0,
         1,1,0,
         0,0,0;

  s_ter<<5,7,5,
         0,0,0,
         0,0,0;

  int sample_size = 10;
  D.resize(sample_size,dimension);

  for(int d=0; d<dimension; d++)
  {
    D.col(d).setLinSpaced(sample_size,s_ini(0,d),s_ter(0,d));
  }

  Eigen::Vector2d s_range = S.get_available_s_range();
  DS.setLinSpaced(sample_size,s_range[0],s_range[1]);
  S.init_with_approximation(s_ini,s_ter,D,DS);

  ceres::Problem problem;

  S.add_vel_cost(&problem);
  S.add_acc_cost(&problem);
  S.add_jerk_cost(&problem);
  S.add_start_cost(&problem);
  S.add_finish_cost(&problem);
  S.add_time_cost(&problem);

  //Set upper and lower boundary
  problem.SetParameterLowerBound(&S.m_beta,0,0.1);
  problem.SetParameterUpperBound(&S.m_beta,0,10.0);
//  problem.SetParameterBlockConstant(&S.m_beta);
  std::cout<<"-------------FINISH PROBLEM CONSTRUCTION------------"<<std::endl;

  //Solve the problem
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout<<summary.FullReport()<<std::endl;
  std::cout<<"-------------FINISH PROBLEM SOLVING------------"<<std::endl;

  //Show the trajectory
  Eigen::Vector2d tr = S.get_available_t_range();
  std::vector<double> t;

  for(double time = tr(0); time <= tr(1); time += 0.1)
  {
    t.push_back(time);
  }

  std::cout<<"------------------"<<std::endl;
  std::cout<<S.get_trajectory(t)<<std::endl;
  std::cout<<"------------------"<<std::endl;
  std::cout<<S.m_beta<<std::endl;

  return 0;
}


//using ceres::AutoDiffCostFunction;
//using ceres::CostFunction;
//using ceres::Problem;
//using ceres::Solve;
//using ceres::Solver;

//struct CostFunctor {
//   template <typename T>
//   bool operator()(const T* const x, T* residual) const {
//     residual[0] = 10.0 - x[0];
//     return true;
//   }
//};

//int main(int argc, char** argv) {
//  google::InitGoogleLogging(argv[0]);

//  // The variable to solve for with its initial value.
//  double initial_x = 5.0;
//  double x = initial_x;

//  // Build the problem.
//  Problem problem;

//  // Set up the only cost function (also known as residual). This uses
//  // auto-differentiation to obtain the derivative (jacobian).
//  CostFunction* cost_function =
//      new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
//  problem.AddResidualBlock(cost_function, nullptr, &x);

//  // Run the solver!
//  Solver::Options options;
//  options.linear_solver_type = ceres::DENSE_QR;
//  options.minimizer_progress_to_stdout = true;
//  Solver::Summary summary;
//  Solve(options, &problem, &summary);

//  std::cout << summary.BriefReport() << "\n";
//  std::cout << "x : " << initial_x
//            << " -> " << x << "\n";
//  return 0;
//}
