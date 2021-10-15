function success = hard_optimization()
global S;

% Fix the initial and terminal conditions
S = S.set_ini_ter_matrix();
S.ctrl_points(1:3,:) = S.A_ini\S.s_ini;
S.ctrl_points(S.n-1:S.n+1,:) = S.A_ter\S.s_ter;

x0 = S.ctrl_points(4:end-3,:);

% The non-linear least square problem
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',100);
% options.Algorithm = 'levenberg-marquardt';

[x,~,~,exit_flag] = lsqnonlin(@mycost,x0,[],[],options);

if (exit_flag > 0)
    S.ctrl_points(4:end-3,:) = reshape(x,S.n+1-6,S.D);
    success = true;
else
    success = false;
end