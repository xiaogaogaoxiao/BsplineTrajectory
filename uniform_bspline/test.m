clear all;
close all;
clc;
N = 40;
beta = 0.20;
M = 30;
S = UniformBspline;
S = S.init(3,N,beta,2);

ctrl_points = 2*rand(N,2);
S = S.set_control_points(ctrl_points);
tr = S.get_available_t_range();
sr = S.get_available_s_range();

% 
trajectory = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],trajectory);
hold on;
s = linspace(tr(1), tr(2), M);
S = S.construct_evaluation_matrix_num(M);
at = S.Eva*S.ctrl_points;
plot(s,at);
% 
nd_traj = (trajectory(2:end,:)-trajectory(1:end-1,:))/0.1;
nd_traj = [nd_traj(1,:); nd_traj];
plot([tr(1):0.1:tr(2)],nd_traj);
% 
dS = S.get_derivative();
% dS = dS.get_derivative();
% % dS = dS.get_derivative();
d_trajectory = dS.get_trajectory([tr(1):0.1:tr(2)]);
% dS = dS.construct_evaluation_matrix_num(20);
% d_at = dS.Eva*dS.ctrl_points;
plot([tr(1):0.1:tr(2)],d_trajectory);
% plot(s,d_at);
% 
S = S.calc_Q_v();
S = S.calc_Q_a();
S = S.calc_Q_j();
% 
% cost = 0;
% for i=1:length(d_trajectory)
%     cost = cost + d_trajectory(i)^2*0.001;
% end
% cost
% 
% % S.beta^5*S.ctrl_points'*S.Q_j*S.ctrl_points
% 

S = S.set_ini_ter_matrix();
S = S.init_with_approximation([1 0 0; 1 0 0]',[5 0 0; 5 0 0]',[1 2 -3 5; 1 5 0 -5]',[8 15 23 31]');
trajectory = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],trajectory);

[F, J] = S.get_vel_cost(1);

% % 
% sx = ([8 15 23 31] - 3)/S.beta;
% sy = [1 2 -3 5];
% plot(sx,sy,'rx');
% % B=[];
% % for s=0:0.1:2
% %     val = S.B_func(1,s);
% %     B=[B;val];
% % end
% % 
% % plot([0:0.1:2],B);