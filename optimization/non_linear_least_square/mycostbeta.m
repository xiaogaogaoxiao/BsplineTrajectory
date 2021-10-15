function [F,J] = mycostbeta(x)
global S;
S.ctrl_points(4:end-3,:) = x(1:end-1);
S.beta = x(end);

[F_smooth, J_smooth] = S.get_jerk_cost_hard(1);

F_time = 10/S.beta;
J_time = -10/S.beta^2;

F = [F_smooth;F_time];
J = [J_smooth;[zeros(1,size(J_smooth,2)-1) J_time]];

J = J(:,[4:end-4, end]);

J(end)
S.beta

