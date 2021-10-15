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
D = 1;
S = UniformBspline;
S = S.init(3,N,beta,D);

S = S.calc_Q_v();
S = S.calc_Q_a();
S = S.calc_Q_j();
S = S.set_ini_ter_matrix();

% The initial and end condition, format is [p; v; a]
s_ini = [0 3 1;  ]';
s_ter = [10 0 0; ]';

% Construct the initial guess
tr = S.get_available_t_range();
sr = S.get_available_s_range();
s = linspace(sr(1),sr(2),10)';

for i=1:D
   d(:,i) =  linspace(s_ini(1,i),s_ter(1,i),10)';
end
S = S.init_with_approximation(s_ini,s_ter,d,s);
init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],init_guess);hold on;

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
figure;
plot([tr(1):0.1:tr(2)],traj);hold on;

dS = S.get_derivative();
traj = dS.get_trajectory([tr(1):0.1:tr(2)+0.1]);
plot([tr(1):0.1:tr(2)+0.1],traj(:,1));hold on;

dS = dS.get_derivative();
traj = dS.get_trajectory([tr(1):0.1:tr(2)+0.1]);
plot([tr(1):0.1:tr(2)+0.1],traj(:,1));hold on;