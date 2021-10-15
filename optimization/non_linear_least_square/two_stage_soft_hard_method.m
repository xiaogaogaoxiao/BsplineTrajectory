clear all;
close all;
clc;

addpath('../../uniform_bspline');

% First step is to initialize the trajectory
global S obs obs2;
obs=[1, 0.9];
obs2=[6.5, 7.5];
N = 20;
beta = 1;
D = 2;
S = UniformBspline;
S = S.init(3,N,beta,D);

% The initial and end condition, format is [p; v; a]
s_ini = [0 0 0;0 0 0]';
s_ter = [10 0 0;10 0 0]';

% Set V,A,J limit
S.v_max = [1, 1, 1]*2;
S.v_min = [-1,-1,-1]*2;

S.a_max = [1, 1, 1]*2;
S.a_min = [-1,-1,-1]*2;

S.j_max = [1, 1, 1]*2;
S.j_min = [-1,-1,-1]*2;

% Construct the initial guess
tr = S.get_available_t_range();
sr = S.get_available_s_range();
s = linspace(sr(1),sr(2),10)';
for i=1:D
   d(:,i) =  linspace(s_ini(1,i),s_ter(1,i),10)';
end
S = S.init_with_approximation(s_ini,s_ter,d,s);

% Plot the initial guess
init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],init_guess);hold on;

% Conduct the soft opmization
soft_optimization();

% Plot the soft optimization's result
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

figure;
plot(traj(:,1),traj(:,2));axis equal;hold on;
plot(init_guess(:,1),init_guess(:,2));axis equal;


circle(obs(1),obs(2),1);
circle(obs2(1),obs2(2),1);
figure;
dS = S.get_derivative();
traj = dS.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

ddS = dS.get_derivative();
traj = ddS.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

dddS = ddS.get_derivative();
traj = dddS.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;


% Conduct the hard opmization
hard_optimization();

% Plot the soft optimization's result
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
figure(2);hold on;
plot(traj(:,1),traj(:,2));axis equal;hold on;