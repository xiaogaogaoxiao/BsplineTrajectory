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

for octt = 1:100000
    % The initial and end condition, format is [p; v; a]
    s_ini = [(rand-0.5)*10 rand-0.5 rand-0.5]'*5;
    s_ter = [(rand-0.5)*20 rand-0.5 rand-0.5]'*-5;
    
    % Construct the initial guess
    tr = S.get_available_t_range();
    sr = S.get_available_s_range();
    s = linspace(sr(1),sr(2),10)';
    d = linspace(s_ini(1),s_ter(1),10)';
    S = S.init_with_approximation(s_ini,s_ter,d,s);
%     init_guess = S.get_trajectory([tr(1):0.1:tr(2)]);
    %     plot([tr(1):0.1:tr(2)],init_guess);hold on;
    x0 = [reshape(S.ctrl_points,[],1); S.beta];
    lb = -inf*ones(size(x0,1),1);
    lb(end)=0.1;
    ub = inf*ones(size(x0,1),1);
    ub(end)=20;

    options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','none');
    % [F,J] = mycost(x0);
    [x,~,~,flag] = lsqnonlin(@mycostbeta_soft,x0,lb,ub,options);
    
    if flag <= 0
        time = datestr(now, 'yyyy_mm_dd');
        filename = sprintf('Optimization_%s.mat',time);
        save(filename);
    end
    octt
end