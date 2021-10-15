% To read more on the open B-spline, please refer to:
% https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-open.html
% and also the paper "Optimal curve fitting and smoothing using normalized
% uniform B-splines: a tool for studying complex systems"
classdef UniformBspline
    properties
        p % Degree
        n % n+1 is the number of control points
        m % m+1 is the number of knots, there is m = n + p +1;
        u % The uniform knot vector [u0, u1, ..., um], u_{i+1} - u_i = 1
        % the valid domain is [up u_{m-p}]
        beta % Time scale t(real time) * beta = u
        ctrl_points % In n+1 x D format
        D % Dimension of the control points
        Eva % Evaluation matrix
        optimize_beta % A flag to set whether to optimize beta, default true
        
        % The following only works for the case p = 3
        Q_v % beta*ctrl_points'*Q_v*ctrl_points is the /int v^2 dt
        Q_a % beta^3*ctrl_points'*Q_v*ctrl_points is the /int a^2 dt
        Q_j % beta^5*ctrl_points'*Q_v*ctrl_points is the /int j^2 dt
        A_ini % initial condition matrix A_ini*[cp1,cp2,cp3]' = [p0 v0 a0]
        A_ter % terminal condition matrix
        s_ini % initial state
        s_ter % terminal state
        
        % Velocity, acceleration and jerk limit
        v_max
        v_min
        a_max
        a_min
        j_max
        j_min
    end
    
    methods
        % Initialization of the uniform b-spline
        function this = init(this,p,n,beta,D)
            this.p = p;
            this.n = n - 1;
            this.beta = beta;
            this.m = this.p + this.n + 1;
            this.u = 0:1:this.m;
            this.D = D;
            this.optimize_beta = true;
            this = this.calc_Q_v();
            this = this.calc_Q_a();
            this = this.calc_Q_j();
        end
        %------------------------------------------------------------------
        % Set the control points
        function this = set_control_points(this, ctrl_points)
            if size(ctrl_points,1) ~= this.n + 1
                error('The size of the control point is wrong!');
            end
            
            if this.D ~= size(ctrl_points,2)
                error('The dimension of the control point is wrong!');
            end
            
            this.ctrl_points = ctrl_points;
            
        end
        %------------------------------------------------------------------
        % Generate the trajectory by evaluting the deboor value on t
        % seconds
        function trajectory = get_trajectory(this, t)
            % t is the vector of time points to evaluate the trajectory
            t_size = length(t);
            
            for t_ctt=1:t_size
                % Map t(t_ctt) to uniform knot vector
                u_probe = t(t_ctt) * this.beta + this.u(this.p + 1);
                trajectory(t_ctt,:) = this.single_deboor(u_probe);
            end
        end
        %------------------------------------------------------------------
        % The deboor's algorithm
        function value = single_deboor(this, u_probe)
            % bound the u_proble
            u_probe = min(max(this.u(this.p + 1), u_probe), this.u(this.m-this.p + 1));
            
            % find which knot segment u_probe belongs to
            k = this.p;
            
            while (true)
                
                if k + 1 + 1 > length(this.u)
                    error('Out side the knot range!')
                end
                
                if this.u(k + 1 + 1) >= u_probe
                    break;
                end
                k = k+1;
            end
            
            % t(t_ctt) is maped to knot segment u_k ~ u_{k+1}
            % there are at most p+1 basis functions N_k-p,p(u), N_k-p+1,p(u), ..., N_k,p(u) non-zero on knot span [uk, uk+1)
            % the effective control points are P_k-p ~ P_k
            % since Matlab index start from 1 instead of 0
            % The effective control points are
            d = this.ctrl_points(k-this.p+1 : k+1, :);
            
            for r=1:this.p
                for i = this.p:-1:r
                    alpha = (u_probe - this.u(i + k - this.p + 1)) / (this.u(i + 1 + k - r + 1) - this.u(i + k - this.p + 1));
                    d(i + 1,:) = (1 - alpha) * d(i - 1 + 1,:) + alpha * d(i + 1,:);
                end
            end
            value = d(this.p + 1,:);
        end
        %------------------------------------------------------------------
        % Construct the evaluation matrix on the provided s range
        function this = construct_evaluation_matrix(this,s)
            this.Eva = zeros(length(s),this.n+1);
            
            for i=1:length(s)
                for j=1:this.n+1
                    this.Eva(i,j) = this.B_func(this.p, s(i)-this.u(j));
                end
            end
            
        end
        %------------------------------------------------------------------
        function range = get_available_s_range(this)
            range = [this.u(this.p+1),this.u(this.m-this.p+1)];
        end
        %------------------------------------------------------------------
        function range = get_available_t_range(this)
            range = [0,this.u(this.m-this.p+1)-this.u(this.p+1)]/this.beta;
        end
        %------------------------------------------------------------------
        function this = construct_evaluation_matrix_num(this,n)
            s = linspace(this.u(this.p+1),this.u(this.m-this.p+1),n);
            this.Eva = zeros(length(s),this.n+1);
            
            for i=1:length(s)
                for j=1:this.n+1
                    this.Eva(i,j) = this.B_func(this.p, s(i)-this.u(j));
                end
            end
            
        end
        %------------------------------------------------------------------
        function val = B_func(this,k,s)
            s_lb = floor(s);
            if s_lb >=0 && s_lb <= k
                val = this.N_func(k-s_lb,k,s-s_lb);
            else
                val = 0;
            end
        end
        %------------------------------------------------------------------
        function pt = N_func(this,j, k, s)
            if j==0 && k == 0
                pt = 1;
            elseif j==0 && k~=0
                pt = (1-s)/k*this.N_func(0,k-1,s);
            elseif j==k && k~=0
                pt = s/k*this.N_func(j-1,k-1,s);
            else
                pt = (k-j+s)/k*this.N_func(j-1,k-1,s) + (1 + j -s)/k*this.N_func(j,k-1,s);
            end
        end
        %------------------------------------------------------------------
        function spline = get_derivative(this)
            spline = this;
            spline.p = spline.p - 1;
            spline.n = spline.n - 1;
            spline.m = spline.p + spline.n + 1;
            spline.u = spline.u(2:end-1);
            
            spline.ctrl_points  = [];
            for i=1:this.n
                spline.ctrl_points(i,:) = this.beta*(this.ctrl_points(i+1,:) - this.ctrl_points(i,:));
            end
        end
        %------------------------------------------------------------------
        function this = calc_Q_v(this)
            if this.p ~= 3
                error('Here we have only pre-calcuated the matrix for 3rd order b-spline, for any other order please use deduction_of_intgration_matrix.mn to calculate it yourself.')
            end
            
            nm = this.n + 1;
            this.Q_v = zeros(nm, nm);
            
            for i=1:nm
                if i == 1
                    this.Q_v(i,i:i+3) = [6 7 -12 -1]/120;
                elseif i == 2
                    this.Q_v(i,i-1:i+3) = [7 40 -22 -24 -1]/120;
                elseif i == 3
                    this.Q_v(i,i-2:i+3) = [-12 -22 74 -15 -24 -1]/120;
                elseif i == nm
                    this.Q_v(i,i-3:i) = [-1 -12 7 6]/120;
                elseif i == nm - 1
                    this.Q_v(i,i-3:i+1) = [-1 -24 -22 40 7]/120;
                elseif i == nm - 2
                    this.Q_v(i,i-3:i+2) = [-1 -24 -15 74 -22 -12]/120;
                else
                    this.Q_v(i,i-3:i+3) = [-1 -24 -15 80 -15 -24 -1]/120;
                end
            end
        end
        %------------------------------------------------------------------
        function this = calc_Q_a(this)
            if this.p ~= 3
                error('Here we have only pre-calcuated the matrix for 3rd order b-spline, for any other order please use deduction_of_intgration_matrix.mn to calculate it yourself.')
            end
            
            nm = this.n + 1;
            this.Q_a = zeros(nm, nm);
            
            for i=1:nm
                if i == 1
                    this.Q_a(i,i:i+3) = [2 -3 0 1]/6;
                elseif i == 2
                    this.Q_a(i,i-1:i+3) = [-3 8 -6 0 1]/6;
                elseif i == 3
                    this.Q_a(i,i-2:i+3) = [0 -6 14 -9 0 1]/6;
                elseif i == nm
                    this.Q_a(i,i-3:i) = [1 0 -3 2]/6;
                elseif i == nm - 1
                    this.Q_a(i,i-3:i+1) = [1 0 -6 8 -3]/6;
                elseif i == nm - 2
                    this.Q_a(i,i-3:i+2) = [1 0 -9 14 -6 0]/6;
                else
                    this.Q_a(i,i-3:i+3) = [1 0 -9 16 -9 0 1]/6;
                end
            end
        end
        %------------------------------------------------------------------
        function this = calc_Q_j(this)
            if this.p ~= 3
                error('Here we have only pre-calcuated the matrix for 3rd order b-spline, for any other order please use deduction_of_intgration_matrix.mn to calculate it yourself.')
            end
            
            nm = this.n + 1;
            this.Q_j = zeros(nm, nm);
            
            for i=1:nm
                if i == 1
                    this.Q_j(i,i:i+3) = [1 -3 3 -1];
                elseif i == 2
                    this.Q_j(i,i-1:i+3) = [-3 10 -12 6 -1];
                elseif i == 3
                    this.Q_j(i,i-2:i+3) = [3 -12 19 -15 6 -1];
                elseif i == nm
                    this.Q_j(i,i-3:i) = [-1 3 -3 1];
                elseif i == nm - 1
                    this.Q_j(i,i-3:i+1) = [-1 6 -12 10 -3];
                elseif i == nm - 2
                    this.Q_j(i,i-3:i+2) = [-1 6 -15 19 -12 3];
                else
                    this.Q_j(i,i-3:i+3) = [-1 6 -15 20 -15 6 -1];
                end
            end
        end
        %------------------------------------------------------------------
        function this = set_ini_ter_matrix(this)
            this.A_ini = [1/6 2/3 1/6;
                [-1/2 0 1/2]*this.beta;
                [1 -2 1]*this.beta^2];
            this.A_ter = [1/6 2/3 1/6;
                [-1/2 0 1/2]*this.beta;
                [1 -2 1]*this.beta^2];
        end
        %------------------------------------------------------------------
        % i is the id of the ctrl point, d is its dimension (both start
        % from 1)
        function id = get_rowid(this,i,d)
            id = (i-1)*this.D + d;
        end
        %------------------------------------------------------------------
        % i is the id of the ctrl point, d is its dimension (both start
        % from 1)
        function id = get_colid(this,i,d)
            id = (d-1)*(this.n+1) + i;
        end
        %------------------------------------------------------------------
        function [F, J] = get_start_cost(this,w)
            b = this.beta;
            M = [1/6 2/3 1/6;
                [-1/2 0 1/2]*b;
                [1 -2 1]*b^2];
            
            F = zeros(3*this.D,1);
            J = zeros(3*this.D,(this.n+1)*this.D + 1);
            
            for d = 1:this.D
                cp = this.ctrl_points(1:3,d);
                idx_a = (d-1)*3 + 1;
                idx_b = d*3;
                F(idx_a:idx_b) = w*(M*cp - this.s_ini(:,d));
                
                J(idx_a,   this.get_colid(1,d)) = w/6;
                J(idx_a+1, this.get_colid(1,d)) = -b*w/2;
                J(idx_a+2, this.get_colid(1,d)) = b^2*w;
                
                J(idx_a,   this.get_colid(2,d)) = 2*w/3;
                J(idx_a+1, this.get_colid(2,d)) = 0;
                J(idx_a+2, this.get_colid(2,d)) = -2*b^2*w;
                
                J(idx_a,   this.get_colid(3,d)) = w/6;
                J(idx_a+1, this.get_colid(3,d)) = b*w/2;
                J(idx_a+2, this.get_colid(3,d)) = b^2*w;
                
                J(idx_a,   end) = 0;
                J(idx_a+1, end) = -w*(cp(1) - cp(3))/2;
                J(idx_a+2, end) = 2*b*w*(cp(1) - 2*cp(2) + cp(3));
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_finish_cost(this,w)
            b = this.beta;
            M = [1/6 2/3 1/6;
                [-1/2 0 1/2]*b;
                [1 -2 1]*b^2];
            
            F = zeros(3*this.D,1);
            J = zeros(3*this.D,(this.n+1)*this.D + 1);
            
            for d = 1:this.D
                cp = this.ctrl_points(end-2:end,d);
                idx_a = (d-1)*3 + 1;
                idx_b = d*3;
                F(idx_a:idx_b) = w*(M*cp - this.s_ter(:,d));
                
                J(idx_a,   this.get_colid(this.n-1,d)) = w/6;
                J(idx_a+1, this.get_colid(this.n-1,d)) = -b*w/2;
                J(idx_a+2, this.get_colid(this.n-1,d)) = b^2*w;
                
                J(idx_a,   this.get_colid(this.n,d)) = 2*w/3;
                J(idx_a+1, this.get_colid(this.n,d)) = 0;
                J(idx_a+2, this.get_colid(this.n,d)) = -2*b^2*w;
                
                J(idx_a,   this.get_colid(this.n+1,d)) = w/6;
                J(idx_a+1, this.get_colid(this.n+1,d)) = b*w/2;
                J(idx_a+2, this.get_colid(this.n+1,d)) = b^2*w;
                
                J(idx_a,   end) = 0;
                J(idx_a+1, end) = -w*(cp(1) - cp(3))/2;
                J(idx_a+2, end) = 2*b*w*(cp(1) - 2*cp(2) + cp(3));
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_obs_cost(this, obs, r, w)
            F = zeros((this.n+1),1);
            J = zeros((this.n+1),(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1
                diff = this.ctrl_points(i,:) - obs;
                dist = norm(diff);
                if (dist < r+0.4)
                    F(i) = (r+0.4-dist)*w;
                    for d=1:this.D
                        J(i, this.get_colid(i,d)) = -w*diff(d)/dist;
                    end
                end
            end
        end
        %------------------------------------------------------------------
        %TODO: It only works for 2D, but since the Matlab version
        %is mainly used for demonstration and code analysis, keep it like
        %this for now.
        function [F, J] = get_map_cost(this, EDT, step,r, w)
            F = zeros((this.n+1),1);
            J = zeros((this.n+1),(this.n+1)*this.D + 1);
            
            % this.n+1 is the total number of control points
            for i=1:this.n+1
                [lb_idx, lb_idy] = this.pos2grid_floor(this.ctrl_points(i,1), this.ctrl_points(i,2),step);
                [prx, pry] = this.pos2grid_raw(this.ctrl_points(i,1), this.ctrl_points(i,2),step);
                
                % bi-linear distance interpolation
                dist = 0;
                for corner_x = lb_idx:lb_idx+1
                    for corner_y = lb_idy:lb_idy+1
                        tmp_val = this.getEDT(EDT,corner_x,corner_y);
                        dist = dist + this.calculate_weight(prx,pry,corner_x,corner_y)*tmp_val;
                    end
                end
                
                if (dist < r+0.05)
                    F(i) = (r+0.05-dist)*w;
                    
                    % Calcuate the gradient by interpolation
                    v00 = this.getEDT(EDT, lb_idx,   lb_idy);
                    v10 = this.getEDT(EDT, lb_idx+1, lb_idy);
                    v11 = this.getEDT(EDT, lb_idx+1, lb_idy+1);
                    v01 = this.getEDT(EDT, lb_idx,   lb_idy+1);
                    % Interpolate along y to get gradient on x
                    k = pry - lb_idy;
                    grad_x = (k*(v11-v01) + (1-k)*(v10-v00))/step;
                    
                    % Interploate along x to get gradient on y
                    k = prx - lb_idx;
                    grad_y = (k*(v11-v10) + (1-k)*(v01-v00))/step;
                    
                    % Assign the gradient
                    J(i, this.get_colid(i,1)) = -w*grad_x;
                    J(i, this.get_colid(i,2)) = -w*grad_y;
                end
            end
        end
        
        function edt = getEDT(this, EDT,idx,idy)
            if idx >=1 && idx <=200 && idy >=1 && idy <=200
                edt = EDT(idx,idy);
            else
                edt = 1000;
            end
        end
        
        function [idx,idy] = pos2grid(this,x,y,step)
            idx = round(x/step)+1;
            idy = round(y/step)+1;
        end
        
        function [idx,idy] = pos2grid_floor(this,x,y,step)
            idx = floor(x/step)+1;
            idy = floor(y/step)+1;
        end
        
        function [idx,idy] = pos2grid_raw(this,x,y,step)
            idx = (x/step)+1;
            idy = (y/step)+1;
        end
        
        function val = calculate_weight(this, prx, pry, corner_x, corner_y)
            val = (1-abs(prx - corner_x))*(1-abs(pry - corner_y));
        end
        %------------------------------------------------------------------
        function [F, J] = get_vel_cost(this, w)
            F = zeros((this.n+1-1)*this.D,1);
            J = zeros((this.n+1-1)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-1
                q = this.ctrl_points(i+1,:) - this.ctrl_points(i,:);
                vel = q*this.beta;
                c = this.beta*w;
                for d=1:this.D
                    F(this.get_rowid(i,d)) = vel(d)*w;
                    J(this.get_rowid(i,d), this.get_colid(i,d)) = -c;
                    J(this.get_rowid(i,d), this.get_colid(i+1,d)) = c;
                    J(this.get_rowid(i,d), end) = q(d)*w;
                end
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_max_vel_cost(this, w)
            F = zeros((this.n+1-1)*this.D,1);
            J = zeros((this.n+1-1)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-1
                q = this.ctrl_points(i+1,:) - this.ctrl_points(i,:);
                vel = q*this.beta;
                c = this.beta*w;
                for d=1:this.D
                    if vel(d) > this.v_max(d)
                        F(this.get_rowid(i,d)) = (vel(d) - this.v_max(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d)) = -c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = c;
                        J(this.get_rowid(i,d), end) = q(d)*w;
                    elseif vel(d) < this.v_min(d)
                        F(this.get_rowid(i,d)) = (-vel(d) + this.v_min(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d)) = c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = -c;
                        J(this.get_rowid(i,d), end) = -q(d)*w;
                    end
                end
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_acc_cost(this, w)
            F = zeros((this.n+1-2)*this.D,1);
            J = zeros((this.n+1-2)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-2
                q = this.ctrl_points(i+2,:) - 2*this.ctrl_points(i+1,:) + this.ctrl_points(i,:);
                acc = q*this.beta^2;
                c = this.beta^2*w;
                for d=1:this.D
                    F(this.get_rowid(i,d)) = acc(d)*w;
                    J(this.get_rowid(i,d), this.get_colid(i,d))   =  c;
                    J(this.get_rowid(i,d), this.get_colid(i+1,d)) = -2*c;
                    J(this.get_rowid(i,d), this.get_colid(i+2,d)) = c;
                    J(this.get_rowid(i,d), end) = 2*w*this.beta*q(d);
                end
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_max_acc_cost(this, w)
            F = zeros((this.n+1-2)*this.D,1);
            J = zeros((this.n+1-2)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-2
                q = this.ctrl_points(i+2,:) - 2*this.ctrl_points(i+1,:) + this.ctrl_points(i,:);
                acc = q*this.beta^2;
                c = this.beta^2*w;
                for d=1:this.D
                    if acc(d) > this.a_max(d)
                        F(this.get_rowid(i,d)) = (acc(d) - this.a_max(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d))   =  c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = -2*c;
                        J(this.get_rowid(i,d), this.get_colid(i+2,d)) = c;
                        J(this.get_rowid(i,d), end) = 2*w*this.beta*q(d);
                    elseif acc(d) < this.a_min(d)
                        F(this.get_rowid(i,d)) = (-acc(d) + this.a_min(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d))   =  -c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = 2*c;
                        J(this.get_rowid(i,d), this.get_colid(i+2,d)) = -c;
                        J(this.get_rowid(i,d), end) = -2*w*this.beta*q(d);
                    end
                end
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_jerk_cost(this, w)
            F = zeros((this.n+1-3)*this.D,1);
            J = zeros((this.n+1-3)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-3
                q = this.ctrl_points(i+3,:) - 3*this.ctrl_points(i+2,:) + 3*this.ctrl_points(i+1,:) - this.ctrl_points(i,:);
                jerk = q*this.beta^3;
                c = this.beta^3*w;
                
                for d=1:this.D
                    F(this.get_rowid(i,d)) = jerk(d)*w;
                    J(this.get_rowid(i,d), this.get_colid(i,d))   =  -c;
                    J(this.get_rowid(i,d), this.get_colid(i+1,d)) = 3*c;
                    J(this.get_rowid(i,d), this.get_colid(i+2,d)) = -3*c;
                    J(this.get_rowid(i,d), this.get_colid(i+3,d)) = c;
                    J(this.get_rowid(i,d), end) = 3*w*this.beta^2*q(d);
                end
            end
        end
        %------------------------------------------------------------------
        function [F, J] = get_max_jerk_cost(this, w)
            F = zeros((this.n+1-3)*this.D,1);
            J = zeros((this.n+1-3)*this.D,(this.n+1)*this.D + 1);
            % this.n+1 is the total number of control points
            for i=1:this.n+1-3
                q = this.ctrl_points(i+3,:) - 3*this.ctrl_points(i+2,:) + 3*this.ctrl_points(i+1,:) - this.ctrl_points(i,:);
                jerk = q*this.beta^3;
                c = this.beta^3*w;
                
                for d=1:this.D
                    if jerk(d) > this.j_max(d)
                        F(this.get_rowid(i,d)) = (jerk(d) - this.j_max(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d))   =  -c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = 3*c;
                        J(this.get_rowid(i,d), this.get_colid(i+2,d)) = -3*c;
                        J(this.get_rowid(i,d), this.get_colid(i+3,d)) = c;
                        J(this.get_rowid(i,d), end) = 3*w*this.beta^2*q(d);
                    elseif jerk(d) < this.j_min(d)
                        F(this.get_rowid(i,d)) = (-jerk(d) + this.j_min(d))*w;
                        J(this.get_rowid(i,d), this.get_colid(i,d))   =  c;
                        J(this.get_rowid(i,d), this.get_colid(i+1,d)) = -3*c;
                        J(this.get_rowid(i,d), this.get_colid(i+2,d)) = 3*c;
                        J(this.get_rowid(i,d), this.get_colid(i+3,d)) = -c;
                        J(this.get_rowid(i,d), end) = -3*w*this.beta^2*q(d);
                    end
                end
            end
        end
        %------------------------------------------------------------------
        % TODO: actually most of the calculation here can be precomputed
        % for increased efficiency.
        function this = init_with_approximation(this,s_ini,s_ter,D,S)
            % Set the init and terminal matrix
            this = set_ini_ter_matrix(this);
            
            % first determine the first and last three control points from the
            % initial and terminal conditions
            this.ctrl_points=zeros(this.n+1,this.D);
            fix = zeros(this.n+1,1);
            %s_ini = this.A_ini * this.ctrl_points(1:3,:);
            this.ctrl_points(1:3,:) = this.A_ini\s_ini;
            fix(1:3)=1;
            
            this.s_ini = s_ini;
            this.s_ter = s_ter;
            
            %similarly for s_ter
            this.ctrl_points(this.n-1:this.n+1,:) = this.A_ter\s_ter;
            fix(this.n-1:this.n+1) = 1;
            
            this = approximation(this,D,S,fix);
        end
        %------------------------------------------------------------------
        % D is the desired data points
        % S is the corrsponding knot point for D
        % fix tells which of the control points shall be fixed
        function this = approximation(this,D,S,fix)
            % Construct the mapping matrix from fix
            k=0;
            M = zeros(this.n+1,this.n+1);
            for i=1:this.n+1
                if fix(i) == 1
                    k = k+1;
                    M(k,i) = 1;
                end
            end
            k_fix = k;
            
            % Get C_f for later use (C_f is the fixed the control points)
            C_f = M*this.ctrl_points;
            C_f = C_f(1:k_fix,:);
            
            % The second scan
            for i=1:this.n+1
                if fix(i) == 0
                    k = k+1;
                    M(k,i) = 1;
                end
            end
            M=M';
            
            %Construct the rest of the matrices
            this = construct_evaluation_matrix(this,S);
            H = this.Eva'*this.Eva + this.beta*this.Q_v*0.5 + this.beta^3*this.Q_a*0.5;
            F = -2*D'*this.Eva;
            H = M'*H*M;
            F = F*M;
            H_fv = H(1:k_fix, k_fix+1:end);
            H_vv = H(k_fix+1:end, k_fix+1:end);
            F_v  = F(:,k_fix+1:end);
            C_v  =  H_vv\(-0.5*(2*C_f'*H_fv+F_v)');
            
            this.ctrl_points = M*[C_f;C_v];
        end
        %------------------------------------------------------------------
    end
end