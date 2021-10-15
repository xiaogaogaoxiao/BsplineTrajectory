function success = soft_optimization()
global S;

% Construct the initial guess for the soft case
x0 = [reshape(S.ctrl_points,[],1); S.beta];
lb = -inf*ones(size(x0,1),1);
lb(end)=1;
ub = inf*ones(size(x0,1),1);
ub(end)=10;

% The non-linear least square problem
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',100);
[x,~,~,exit_flag] = lsqnonlin(@mycostbeta_soft,x0,lb,ub,options);

if (exit_flag > 0)
    S.ctrl_points = reshape(x(1:end-1),[],S.D);
    S.beta = x(end);
    success = true;
else
    success = false;
end