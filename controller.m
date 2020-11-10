function [u,J_OL] = controller(sys,cont,xk,ref)

%% setup optimization problem

% Declare independent variables
% l goes from 1:N
% j goes from 1:nx_v
z_lk = sdpvar(sys.n,cont.N+1,'full');
v_lk = sdpvar(sys.m,cont.N,'full');
alpha_lk = sdpvar(1,cont.N+1,'full');
Lambda_jlk = sdpvar(cont.nHx,sys.nHtheta,cont.nx_v,cont.N-1,'full');
Lambda_1k = sdpvar(cont.nHx,sys.nHtheta,'full');

% define state dependent variables
x_jlk = [];
u_jlk = [];
d_jlk = [];
D_jlk = [];
for l = 2:cont.N    
    x_jlk = cat(3,x_jlk,repmat(z_lk(:,l),1,cont.nx_v)+alpha_lk(:,l)*cont.x_v);
    u_jlk = cat(3,u_jlk,repmat(v_lk(:,l)-cont.K*ref(l).x_s,1,cont.nx_v)+cont.K*x_jlk(:,:,l-1));
    
    d_jlk = cat(3,d_jlk,sys.A0*x_jlk(:,:,l-1)+sys.B0*u_jlk(:,:,l-1) - repmat(z_lk(:,l+1),1,cont.nx_v));
    D_temp = [];
    for j = 1:cont.nx_v
        D_temp = cat(3,D_temp,D_mult(sys,x_jlk(:,j,l-1),u_jlk(:,j,l-1)));
    end    
    D_jlk = cat(4,D_jlk,D_temp);    
end

%%
% define constraints.
Constraints = [Lambda_jlk(:)>=0,Lambda_1k(:)>=0];
Constraints = [Constraints,z_lk(:,1)==xk,alpha_lk(1)==0];
Constraints = [Constraints,...
                (sys.F+sys.G*cont.K)*xk + sys.G*(v_lk(:,1)-cont.K*ref(1).x_s) <= ones(sys.nc,1),...
                Lambda_1k*cont.h_theta_k + cont.H_x*(sys.A0*xk+sys.B0*(cont.K*(xk-ref(1).x_s)+v_lk(:,1))-z_lk(:,2)) - alpha_lk(2)*ones(cont.nHx,1) <= -cont.w_bar,...
                cont.H_x*D_mult(sys,xk,cont.K*(xk-ref(1).x_s)+v_lk(:,1)) == Lambda_1k*cont.H_theta
                ];
for l = 2:cont.N
     Constraints = [Constraints, ...s
         (sys.F+sys.G*cont.K)*z_lk(:,l) + sys.G*(v_lk(:,l)-cont.K*ref(l).x_s) + alpha_lk(:,l)*cont.f_bar <= ones(sys.nc,1)];         
     
     for j = 1:cont.nx_v
         Constraints = [Constraints, ...
             Lambda_jlk(:,:,j,l-1)*cont.h_theta_k + cont.H_x*d_jlk(:,j,l-1) - alpha_lk(:,l+1)*ones(cont.nHx,1) <= -cont.w_bar];
         
         Constraints = [Constraints, ...
            cont.H_x*D_jlk(:,:,j,l-1) == Lambda_jlk(:,:,j,l-1)*cont.H_theta];         
     end
end

% define terminal constraints
Constraints = [Constraints, ...
               z_lk(:,cont.N+1) == ref(end).x_s, ...
                alpha_lk(cont.N+1) <= ref(end).alpha_T...
                alpha_lk(cont.N+1) >= cont.alpha_min];    
            
options = sdpsettings('solver','gurobi','verbose',0,'debug',0);


%% Cost function definition

% define cost variables
stage_cost_max = sdpvar(cont.N+1,1,'full');
J = sum(stage_cost_max);

%% define cost constraints

% define cost variables
stage_cost_max = sdpvar(cont.N+1,1,'full');
J = sum(stage_cost_max);

costConstraints = [norm(cont.R_L*(cont.K*(xk-ref(1).x_s)+v_lk(:,1)),'inf')<=stage_cost_max(1)];

for l = 2:cont.N
    for j = 1:cont.nx_v                
        costConstraints = [costConstraints, ... 
            norm(cont.Q_L*(x_jlk(:,j,l-1)-ref(l).x_s),'inf') +  norm(cont.R_L*u_jlk(:,j,l-1),'inf') <= stage_cost_max(l)
        ];
    end
end

% terminal cost
for j = 1:cont.nx_v
    costConstraints = [costConstraints, ... 
        alpha_lk(cont.N+1)*(norm(cont.Q_L*cont.x_v(:,j),'inf') + norm(cont.R_L*(cont.K*cont.x_v(:,j)+ref(end).u_s),'inf')) <= stage_cost_max(cont.N+1)
    ];
end

diagnostics = optimize([Constraints,costConstraints],J,options);
if diagnostics.problem
   error(diagnostics.info);
end
u = value(cont.K*(xk-ref(1).x_s) + v_lk(:,1));
J_OL = value(J);
end

