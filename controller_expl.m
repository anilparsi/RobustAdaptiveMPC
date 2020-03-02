function u = controller_expl(sys,cont,xk)

%% Part for standard constraints
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
    x_jlk = cat(3,x_jlk,repmat(z_lk(:,l),1,cont.nx_v)+alpha_lk(l)*cont.x_v);
    u_jlk = cat(3,u_jlk,repmat(v_lk(:,l),1,cont.nx_v)+cont.K*x_jlk(:,:,l-1));
    
    d_jlk = cat(3,d_jlk,sys.A0*x_jlk(:,:,l-1)+sys.B0*u_jlk(:,:,l-1) - repmat(z_lk(:,l+1),1,cont.nx_v));
    D_temp = [];
    for j = 1:cont.nx_v
        D_temp = cat(3,D_temp,D_mult(sys,x_jlk(:,j,l-1),u_jlk(:,j,l-1)));
    end    
    D_jlk = cat(4,D_jlk,D_temp);    
end

% define constraints.
stdConstraints = [Lambda_jlk(:)>=0,Lambda_1k(:)>=0];
stdConstraints = [stdConstraints,z_lk(:,1)==xk,alpha_lk(1)==0];
stdConstraints = [stdConstraints,...
                    (sys.F+sys.G*cont.K)*xk + sys.G*v_lk(:,1) <= ones(sys.nc,1),...
                    Lambda_1k*cont.h_theta_k + cont.H_x*(sys.A0*xk+sys.B0*(cont.K*xk+v_lk(:,1))-z_lk(:,2)) - alpha_lk(2)*ones(cont.nHx,1) <= -cont.w_bar,...
                    cont.H_x*D_mult(sys,xk,cont.K*xk+v_lk(:,1)) == Lambda_1k*cont.H_theta
                    ];
for l = 2:cont.N
     stdConstraints = [stdConstraints, ...
         (sys.F+sys.G*cont.K)*z_lk(:,l) + sys.G*v_lk(:,l) + alpha_lk(l)*cont.f_bar <= ones(sys.nc,1)];         
     
     for j = 1:cont.nx_v
         stdConstraints = [stdConstraints, ...
            Lambda_jlk(:,:,j,l-1)*cont.h_theta_k + cont.H_x*d_jlk(:,j,l-1) - alpha_lk(l+1)*ones(cont.nHx,1) <= -cont.w_bar];
         
         stdConstraints = [stdConstraints, ...
            cont.H_x*D_jlk(:,:,j,l-1) == Lambda_jlk(:,:,j,l-1)*cont.H_theta];         
     end
end

% define terminal constraints
stdConstraints = [stdConstraints, ...
               z_lk(:,cont.N+1) == zeros(sys.n,1), ...
                cont.h_T*alpha_lk(cont.N+1) <= 1 ];    

%% Part for exploration

% Declare variables for exploration
if cont.nPred_X>1
    Lambda_til_jlk = sdpvar(cont.nHx,sys.nHtheta+cont.nPred_theta*sys.nHw,cont.nx_v,cont.nPred_X-1,'full');
else
    Lambda_til_jlk = [];
end
Lambda_til_1k = sdpvar(cont.nHx,sys.nHtheta+cont.nPred_theta*sys.nHw);
z_til_lk = sdpvar(sys.n,cont.nPred_X+1,'full');
alpha_til_lk = sdpvar(1,cont.nPred_X+1,'full');

% predict state evolution
x_hat = xk;
for l = 1:cont.nPred_theta-1
    x_hat = [x_hat  (cont.A_est+cont.B_est*cont.K)*x_hat(:,l) + cont.B_est*v_lk(:,l)];
end  

% predict new constraints on theta
H_pred = [sys.H_w*D_mult(sys,xk,cont.K*xk+v_lk(:,1))];
h_pred = [sys.h_w + sys.H_w*D_mult(sys,xk,cont.K*xk+v_lk(:,1))*cont.theta_hat];
for l = 2:cont.nPred_theta
    H_pred = [H_pred; 
              sys.H_w*D_mult(sys,x_hat(:,l),cont.K*x_hat(:,l)+v_lk(:,l))];
    h_pred = [h_pred; 
              sys.h_w + sys.H_w*D_mult(sys,x_hat(:,l),cont.K*x_hat(:,l)+v_lk(:,l))*cont.theta_hat];
end

% Append to H_theta
H_theta_til = [cont.H_theta; H_pred];
h_theta_til = [cont.h_theta_k; h_pred];

% define exploration variables
x_til_jlk = [];
u_til_jlk = [];
d_til_jlk = [];
D_til_jlk = [];

for l = 2:cont.nPred_X+1
    x_til_jlk = cat(3,x_til_jlk,repmat(z_til_lk(:,l),1,cont.nx_v)+alpha_til_lk(l)*cont.x_v);
    u_til_jlk = cat(3,u_til_jlk,repmat(v_lk(:,l),1,cont.nx_v)+cont.K*x_til_jlk(:,:,l-1));
end
for l = 2:cont.nPred_X    
    d_til_jlk = cat(3,d_til_jlk,sys.A0*x_til_jlk(:,:,l-1)+sys.B0*u_til_jlk(:,:,l-1) - repmat(z_til_lk(:,l+1),1,cont.nx_v));
    D_temp = [];
    for j = 1:cont.nx_v
        D_temp = cat(3,D_temp,D_mult(sys,x_til_jlk(:,j,l-1),u_til_jlk(:,j,l-1)));
    end    
    D_til_jlk = cat(4,D_til_jlk,D_temp);    
end

% define exploration constraints
exploreConstraints = [Lambda_til_jlk(:)>=0,Lambda_til_1k(:)>=0];
exploreConstraints = [exploreConstraints,z_til_lk(:,1)==xk,alpha_til_lk(1)==0];
exploreConstraints = [exploreConstraints,...
                    Lambda_til_1k*h_theta_til + cont.H_x*(sys.A0*xk+sys.B0*(cont.K*xk+v_lk(:,1))-z_til_lk(:,2)) - alpha_til_lk(2)*ones(cont.nHx,1) <= -cont.w_bar,...
                    cont.H_x*D_mult(sys,xk,cont.K*xk+v_lk(:,1)) == Lambda_til_1k*H_theta_til
                    ];
% implement containment conditions of tubes
for l = 2:cont.nPred_X
    exploreConstraints = [exploreConstraints, ...
    (sys.F+sys.G*cont.K)*z_til_lk(:,l) + sys.G*v_lk(:,l) + alpha_til_lk(l)*cont.f_bar <= ones(sys.nc,1)
    ];
     for j = 1:cont.nx_v
         exploreConstraints = [exploreConstraints, ...
             Lambda_til_jlk(:,:,j,l-1)*h_theta_til + cont.H_x*d_til_jlk(:,j,l-1) - alpha_til_lk(l+1)*ones(cont.nHx,1) <= -cont.w_bar];
         
         exploreConstraints = [exploreConstraints, ...
            cont.H_x*D_til_jlk(:,:,j,l-1) == Lambda_til_jlk(:,:,j,l-1)*H_theta_til];         
     end
end
% at the end, ensure the X_til_lk tube is inside the X_lk tube
for j = 1:cont.nx_v
    exploreConstraints = [exploreConstraints, ...
        cont.H_x*(z_til_lk(:,end)+alpha_til_lk(end)*cont.x_v(:,j) - z_lk(:,cont.nPred_X+1))  <= alpha_lk(cont.nPred_X+1)*cont.vec_1_x;    ];
end
%% Cost function definition

% define cost variables
stage_cost_max = sdpvar(cont.N+1,1,'full');
J = sum(stage_cost_max);

% define cost constraints
costConstraints = [norm(cont.R_L*(cont.K*xk+v_lk(:,1)),'inf')<=stage_cost_max(1)];

for l = 2:cont.N
    for j = 1:cont.nx_v                
        if l<=cont.nPred_X+1
            costConstraints = [costConstraints, ... 
                norm(cont.Q_L*x_til_jlk(:,j,l-1),'inf') +  norm(cont.R_L*(cont.K*x_til_jlk(:,j,l-1)+v_lk(:,l)),'inf') <= stage_cost_max(l)
            ];
        else
            costConstraints = [costConstraints, ... 
                norm(cont.Q_L*x_jlk(:,j,l-1),'inf') +  norm(cont.R_L*(cont.K*x_jlk(:,j,l-1)+v_lk(:,l)),'inf') <= stage_cost_max(l)
            ];
        end
    end
end
% terminal cost
for j = 1:cont.nx_v
    costConstraints = [costConstraints, ... 
        norm(cont.Q_L*cont.x_v(:,j),'inf') + norm(cont.R_L*cont.K*cont.x_v(:,j),'inf') <= stage_cost_max(cont.N+1)
    ];
end

%% Solve problem
Constraints = [stdConstraints;exploreConstraints;costConstraints];    
options = sdpsettings('solver','ipopt','verbose',0);
diagnostics = optimize(Constraints,J,options);
if diagnostics.problem
   error(diagnostics.info) 
end
u = value(cont.K*xk + v_lk(:,1));

end

