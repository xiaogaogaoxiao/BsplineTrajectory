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
s_ini = [0 0 0]';
s_ter = [5 0 0]';

% Construct the initial guess
tr = S.get_available_t_range();
sr = S.get_available_s_range();
s = linspace(sr(1),sr(2),10)';
d = linspace(s_ini(1),-s_ter(1),10)';
S = S.init_with_approximation(s_ini,s_ter,d,s);
init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],init_guess);hold on;
x0 = S.ctrl_points(4:end-3,:);

% The non-linear least square problem
% Turn on the gradient first
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
options.Algorithm = 'levenberg-marquardt';
% [F,J] = mycost(x0);
[x] = lsqnonlin(@mycost,x0,[],[],options);
% 
x = reshape(x,S.n+1-6,S.D);
S.ctrl_points(4:end-3,:) = x;
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;
