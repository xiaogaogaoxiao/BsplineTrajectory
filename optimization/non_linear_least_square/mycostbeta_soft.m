function [F,J] = mycostbeta_soft(x)
global S EDT; 
S.ctrl_points = reshape(x(1:end-1),[],S.D);
S.beta = x(end);

[Fsv, Jsv] = S.get_vel_cost(1);
[Fsa, Jsa] = S.get_acc_cost(1);
[Fsj, Jsj] = S.get_jerk_cost(1);

[Fmv, Jmv] = S.get_max_vel_cost(5);
[Fma, Jma] = S.get_max_acc_cost(5);
[Fmj, Jmj] = S.get_max_jerk_cost(5);

[F_start, J_start] = S.get_start_cost(32);
[F_end, J_end] = S.get_finish_cost(32);
%--------------------------------------------------------------------------
r = 1.0;
[Fmap, Jmap] = S.get_map_cost(EDT,0.05,r,50);
%--------------------------------------------------------------------------
F_time = 20/S.beta;
J_time = -20/S.beta^2;

F = [Fsv;Fsa;Fsj;Fmv;Fma;Fmj;F_start;F_end;Fmap;F_time];
J = [Jsv;Jsa;Jsj;Jmv;Jma;Jmj;J_start;J_end;Jmap;[zeros(1,size(Jsv,2)-1) J_time]];

% J = J(:,[4:end-4, end]);

