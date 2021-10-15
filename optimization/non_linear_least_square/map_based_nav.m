clear all;
close all;
clc;

global S EDT obs;
% setup the obstacle map
bw = zeros(200,200);
step = 0.05;
obs =[5.1 3;4.9 7;4.5 5.1];
for i=1:size(obs,1)
    [idx,idy] = pos2grid(obs(i,1),obs(i,2),step);
    bw(idx,idy) = 1;
end
EDT = bwdist(bw,'euclidean')*step;


DX = [0:1:size(EDT,1)-1]*step;
DY = [0:1:size(EDT,2)-1]*step;
contourf(DX,DY,EDT');hold on;

addpath('../../uniform_bspline');

% First step is to initialize the trajectory

N = 40;
beta = 4;
D = 2;
S = UniformBspline;
S = S.init(3,N,beta,D);

% The initial and end condition, format is [p; v; a]
s_ini = [5 0 0;1 0 0]';
s_ter = [5 0 0;9 0 0]';

% Set V,A,J limit
S.v_max = [1, 1, 1]*1;
S.v_min = [-1,-1,-1]*1;

S.a_max = [1, 1, 1]*1;
S.a_min = [-1,-1,-1]*1;

S.j_max = [1, 1, 1]*1;
S.j_min = [-1,-1,-1]*1;

% Construct the initial guess
tr = S.get_available_t_range();
sr = S.get_available_s_range();
s = linspace(sr(1),sr(2),10)';
for i=1:D
   d(:,i) =  linspace(s_ini(1,i),s_ter(1,i),10)';
end
S = S.init_with_approximation(s_ini,s_ter,d,s);

% Plot the initial guess
figure(2);
init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],init_guess);hold on;

% Conduct the soft opmization
soft_optimization();

% Plot the soft optimization's result
figure(2);
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

figure(1);
plot(traj(:,1),traj(:,2),'y');axis equal;hold on;
plot(init_guess(:,1),init_guess(:,2));axis equal;

for i=1:size(obs,1)
    circle(obs(i,1),obs(i,2),1);
end


% Conduct the hard opmization
hard_optimization();

% Plot the soft optimization's result
figure(2);
tr = S.get_available_t_range();
traj = S.get_trajectory([tr(1):0.1:tr(2)]);
plot([tr(1):0.1:tr(2)],traj);hold on;

figure(1);hold on;
plot(traj(:,1),traj(:,2));axis equal;hold on;
plot(init_guess(:,1),init_guess(:,2));axis equal;

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

