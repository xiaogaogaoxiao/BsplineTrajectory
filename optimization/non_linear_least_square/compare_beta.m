clear all;
close all;
clc;

addpath('../../uniform_bspline');

% Our first simple test problem is to use the non-linear least square solver
%problem to generate a smooth trajectory from one fixed start condition to
%a fixed end condition with fixed total time.

% First step is to initialize the trajectory
global S;
N = 50;
beta = 10;
S = UniformBspline;
S = S.init(3,N,beta,1);

S = S.calc_Q_v();
S = S.calc_Q_a();
S = S.calc_Q_j();
S = S.set_ini_ter_matrix();

% The initial and end condition, format is [p; v; a]
s_ini = [rand rand rand]'*5;
s_ter = [-100 rand rand]'*-5;

% Construct the initial guess
tr = S.get_available_t_range();
sr = S.get_available_s_range();
s = linspace(sr(1),sr(2),10)';
d = linspace(s_ini(1),s_ter(1),10)';
S = S.init_with_approximation(s_ini,s_ter,d,s);
init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],init_guess);hold on;
x0 = [S.ctrl_points(4:end-3,:); S.beta];
lb = -inf*ones(size(x0,1),1);
lb(end)=0.1;
ub = inf*ones(size(x0,1),1);
ub(end)=20;
% The non-linear least square problem
% Turn on the gradient first
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','iter');
% [F,J] = mycost(x0);
[x] = lsqnonlin(@mycostbeta,x0,lb,ub,options);

S.ctrl_points(4:end-3,:) = x(1:end-1);
S.beta = x(end);
S = S.set_ini_helper();
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.01:tr(2)+0.1]);
plot([tr(1):0.01:tr(2)+0.1],traj);hold on;

%% Here enters the comparison phase
x0 = S.ctrl_points(4:end-3,:);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
options.Algorithm = 'levenberg-marquardt';
[x] = lsqnonlin(@mycost,x0,[],[],options);
x = reshape(x,S.n+1-6,S.D);
S.ctrl_points(4:end-3,:) = x;
traj = S.get_trajectory([tr(1):0.01:tr(2)+0.1]);
plot([tr(1):0.01:tr(2)+0.1],traj);hold on;

%% Compare with soft constraint case
x0 = [reshape(S.ctrl_points,[],1); S.beta];
lb = -inf*ones(size(x0,1),1);
lb(end)=0.1;
ub = inf*ones(size(x0,1),1);
ub(end)=20;
% The non-linear least square problem
% Turn on the gradient first
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','iter');
% [F,J] = mycost(x0);
[x] = lsqnonlin(@mycostbeta_soft,x0,lb,ub,options);

S.ctrl_points = reshape(x(1:end-1),[],S.D);
S.beta = x(end);
% S = S.set_ini_helper();
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

% %% Compare the F values
% [F1, J1] = S.get_acc_cost(1);
% [F2, J2] = S.get_acc_cost_hard(1);
% max(abs(F1-F2))
% %% Compare the F gradient
% dt = 0.001;
% [Fa, Ja] = S.get_jerk_cost_hard(1);
% S.beta = S.beta + dt;
% [Fb, Jb] = S.get_jerk_cost_hard(1);
% 
% J_num = (Fb - Fa)/dt;
% J_ana = Ja(:,end);
% max(abs(J_num-J_ana))
% 
% %%
% [F_start, J_start] = S.get_start_cost(20);
