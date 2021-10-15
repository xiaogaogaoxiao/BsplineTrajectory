function [F,J] = mycost(x)
global S EDT;
S.ctrl_points(4:end-3,:) = reshape(x,[],S.D);

[Fsv, Jsv] = S.get_vel_cost(1);
[Fsa, Jsa] = S.get_acc_cost(1);
[Fsj, Jsj] = S.get_jerk_cost(1);

[Fmv, Jmv] = S.get_max_vel_cost(5);
[Fma, Jma] = S.get_max_acc_cost(5);
[Fmj, Jmj] = S.get_max_jerk_cost(5);

%--------------------------------------------------------------------------
r = 1.0;
[Fmap, Jmap] = S.get_map_cost(EDT,0.05,r,50);
%--------------------------------------------------------------------------

F = [Fsv;Fsa;Fsj;Fmv;Fma;Fmj;Fmap;];
J = [Jsv;Jsa;Jsj;Jmv;Jma;Jmj;Jmap;];

N = S.n+1;
M = zeros(1,(S.n+1-6)*S.D);
for d=0:S.D-1
    shift = d*(S.n+1-6);
    M(shift+1 : shift+S.n+1-6) = [d*N+4 : d*N+S.n-2];
end
 
J = J(:,M);
