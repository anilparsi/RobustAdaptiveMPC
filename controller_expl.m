function u = controller_expl(sys,cont,xk)

%% Part for standard constraints
% Declare independent variables
% l goes from 1:N
% j goes from 1:nx_v
z_lk = sdpvar(sys.n,cont.N+1,'full');
u_lk = sdpvar(sys.m,cont.N,'full');
alpha_lk = sdpvar(1,cont.N+1,'full');
Lambda_jlk = sdpvar(cont.nHx,sys.nHtheta,cont.nx_v,cont.N,'full');

% define state dependent variables
x_jlk = [];
u_jlk = [];
d_jlk = [];
D_jlk = [];
for l = 1:cont.N    
    x_jlk = cat(3,x_jlk,repmat(z_lk(:,l),1,cont.nx_v)+alpha_lk(:,l)*cont.x_v);
    u_jlk = cat(3,u_jlk,repmat(u_lk(:,l),1,cont.nx_v));
    
    d_jlk = cat(3,d_jlk,sys.A0*x_jlk(:,:,l)+sys.B0*u_jlk(:,:,l) - repmat(z_lk(:,l+1),1,cont.nx_v));
    D_temp = [];
    for j = 1:cont.nx_v
        D_temp = cat(3,D_temp,D_mult(sys,x_jlk(:,j,l),u_jlk(:,j,l)));
    end    
    D_jlk = cat(4,D_jlk,D_temp);    
end

% define constraints.
Constraints = [Lambda_jlk(:)>=0];
Constraints = [Constraints,-cont.H_x*z_lk(:,1) - alpha_lk(:,1)*ones(cont.nHx,1) <= -cont.H_x*xk];
for l = 1:cont.N
     Constraints = [Constraints, ...
         sys.F*z_lk(:,l) + sys.G*u_lk(:,l) + alpha_lk(:,l)*cont.f_bar <= ones(sys.nc,1)];         
     
     for j = 1:cont.nx_v
         Constraints = [Constraints, ...
             Lambda_jlk(:,:,j,l)*cont.h_theta_k + cont.H_x*d_jlk(:,j,l) - alpha_lk(:,l+1)*ones(cont.nHx,1) <= -cont.w_bar];
         
         Constraints = [Constraints, ...
            cont.H_x*D_jlk(:,:,j,l) == Lambda_jlk(:,:,j,l)*cont.H_theta];         
     end
end

% define terminal constraints
Constraints = [Constraints, ...
               z_lk(:,cont.N+1) == zeros(sys.n,1), ...
                cont.h_T*alpha_lk(cont.N+1) <= 1 ];    

%% Part for exploration

% Declare variables for exploration
Lambda_til_jlk = sdpvar(cont.nHx,sys.nHtheta+cont.nPred_theta*sys.nHw,cont.nx_v,cont.nPred_X,'full');
z_til_lk = sdpvar(sys.n,cont.nPred_X+1,'full');
alpha_til_lk = sdpvar(1,cont.nPred_X+1,'full');

% predict state evolution
x_hat = xk;
for l = 1:cont.nPred_theta-1
    x_hat = [x_hat  cont.A_est*x_hat(:,l) + cont.B_est*u_lk(:,l)];
end  

% predict new constraints on theta
H_pred = [sys.H_w*D_mult(sys,xk,u_lk(:,1))];
h_pred = [sys.h_w + sys.H_w*D_mult(sys,xk,u_lk(:,1))*cont.theta_hat];
for l = 2:cont.nPred_theta-1
    H_pred = [H_pred; 
              sys.H_w*D_mult(sys,x_hat(:,l),u_hat(:,l))];
    h_pred = [h_pre; 
              sys.h_w + sys.H_w*D_mult(sys,x_hat(:,l),u_hat(:,l))*cont.theta_hat];
end

% Append to H_theta
H_theta_til = [cont.H_theta; H_pred];
h_theta_til = [cont.h_theta_k; h_pred];

% define exploration variables
x_til_jlk = [];
u_til_jlk = [];
d_til_jlk = [];
D_til_jlk = [];
for l = 1:cont.nPred_X   
    x_til_jlk = cat(3,x_til_jlk,repmat(z_til_lk(:,l),1,cont.nx_v)+alpha_til_lk(:,l)*cont.x_v);
    u_til_jlk = cat(3,u_til_jlk,repmat(u_lk(:,l),1,cont.nx_v));
    
    d_til_jlk = cat(3,d_til_jlk,sys.A0*x_til_jlk(:,:,l)+sys.B0*u_til_jlk(:,:,l) - repmat(z_til_lk(:,l+1),1,cont.nx_v));
    D_temp = [];
    for j = 1:cont.nx_v
        D_temp = cat(3,D_temp,D_mult(sys,x_til_jlk(:,j,l),u_til_jlk(:,j,l)));
    end    
    D_til_jlk = cat(4,D_til_jlk,D_temp);    
end

% define exploration constraints
for l = 1:cont.nPred_X
    Constraints = [Constraints, ...
    sys.F*z_til_lk(:,l) + sys.G*u_lk(:,l) + alpha_til_lk(:,l)*cont.f_bar <= ones(sys.nc,1)
    ];
     for j = 1:cont.nx_v
         Constraints = [Constraints, ...
             Lambda_til_jlk(:,:,j,l)*h_theta_til + cont.H_x*d_til_jlk(:,j,l) - alpha_til_lk(:,l+1)*ones(cont.nHx,1) <= -cont.w_bar];
         
         Constraints = [Constraints, ...
            cont.H_x*D_til_jlk(:,:,j,l) == Lambda_til_jlk(:,:,j,l)*H_theta_til];         
     end
end
for j = 1:cont.nx_v
    Constraints = [Constraints, ...
        cont.H_x*(z_til_lk(:,end)+alpha_til_lk(:,end)*cont.x_v(:,j) - z_lk(:,cont.nPred_X+1))  <= alpha_lk(:,cont.nPred_X+1)*cont.vec_1_x;    ];
end
%% Cost function definition

% define cost variables
stage_cost_max_x = sdpvar(sys.n,cont.N+1,'full');
stage_cost_u = sdpvar(sys.m,cont.N,'full');
J = sum(stage_cost_max_x(:,1));
for l = 1:cont.N
    J = J + sum(stage_cost_max_x(:,l+1))+ sum(stage_cost_u(:,l));
end

% define cost constraints
for l = 1:cont.N
    for j = 1:cont.nx_v                
        if l<=cont.nPred_X
            Constraints = [Constraints, ... 
                cont.Q*(z_til_lk(:,l) + alpha_til_lk(l)*cont.x_v(:,j)) <= stage_cost_max_x(:,l), ...
                -cont.Q*(z_til_lk(:,l) + alpha_til_lk(l)*cont.x_v(:,j)) <= stage_cost_max_x(:,l), ...
                cont.R*u_jlk(:,j,l) <= stage_cost_u(:,l),...
                -cont.R*u_jlk(:,j,l) <= stage_cost_u(:,l)
            ];
        else
            Constraints = [Constraints, ... 
                cont.Q*x_jlk(:,j,l) <= stage_cost_max_x(:,l), ...
                -cont.Q*x_jlk(:,j,l) <= stage_cost_max_x(:,l), ...
                cont.R*u_jlk(:,j,l) <= stage_cost_u(:,l),...
                -cont.R*u_jlk(:,j,l) <= stage_cost_u(:,l)
            ];
        end
    end
end
% terminal cost
for j = 1:cont.nx_v
    Constraints = [Constraints, ... 
        (cont.P)*alpha_lk(cont.N+1)*cont.x_v(:,j) <= stage_cost_max_x(:,cont.N+1)
    ];
end

%% Solve problem
    
options = sdpsettings('solver','ipopt');%,'fmincon','fmincon.Maxiter',20
diagnostics = optimize(Constraints,J,options);
if diagnostics.problem
   error(diagnostics.info) 
end
u = value(u_lk(:,1));



end

