function [u,J_OL,warn_flag] = controller_PE(sys,cont,xk,U_past,PE,ref)

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

Constraints = [Constraints;costConstraints];
%% Define persistency of excitation conditions

% Since PE is convex only for single inputs, we choose which input to apply
% through the variable PE
u_hat = cont.K*(xk-ref(1).x_s)+ v_lk(:,1);
omega_bar_PE = -cont.rho_PE*eye(sys.n-1);
for k =0:cont.P_PE-1
    phi_vec = vec(U_past(PE,k+1:k+sys.n-1));
    omega_bar_PE = omega_bar_PE + phi_vec*phi_vec';
end   
phi_vec0 = vec(U_past(PE,1:sys.n-1));


alpha_PE = 1-phi_vec0'*(omega_bar_PE\phi_vec0);
beta_PE = 0;
u_phi_vec = 0;
u_u = 0;
for j = 1:cont.P_PE-1
    phi_vec = vec(U_past(PE,j+1:j+sys.n-1));
    beta_PE = beta_PE - U_past(PE,j)*phi_vec'*(omega_bar_PE\phi_vec0);
    u_phi_vec = u_phi_vec + U_past(PE,j)*phi_vec';
    u_u = u_u + U_past(PE,j)*U_past(PE,j);
end
gamma_PE = u_u - cont.rho_PE...
            - u_phi_vec'*(omega_bar_PE\u_phi_vec);

warn_flag = 0;       
try
    bounds = roots([alpha_PE,2*beta_PE,gamma_PE]);
    bounds = sort(bounds);
catch
    % Past inputs are not PE. Reduce the PE order to 1, use a sufficient
    % condition:
    warn_flag = 1;       
    alpha_PE = 1; % choose a positive number
    bounds = [-sqrt(cont.rho_PE); sqrt(cont.rho_PE)];    
end
         
try
    if alpha_PE < 0
        % stay within the bounds
        Constraints0 = [Constraints,...
                        bounds(1) <= u_hat(PE) <= bounds(2)   ];
        diagnostics = optimize(Constraints0,J,options);
        if diagnostics.problem
           error('Infeasibility due to PE') 
        end
        u = value(u_hat(:,1));
    elseif alpha_PE == 0
        Constraints0 = [Constraints,...
                        2*beta_PE*u_hat(PE) + gamma_PE >=0   ];
        diagnostics = optimize(Constraints0,J,options);
        if diagnostics.problem
           error('Infeasibility due to PE') 
        end
        u = value(u_hat(:,1));    
            
    elseif isreal(bounds(1))
        % stay out of the bounds
        Constraints1 = [Constraints,...
                        u_hat(PE) <= bounds(1)   ];
        Constraints2 = [Constraints,...  
                        u_hat(PE) >= bounds(2) ];


        diagnostics1 = optimize(Constraints1,J,options);
        if ~diagnostics1.problem
            J1 = value(J);
            u1 = value(u_hat(:,1));
        else
            J1 = Inf;
            u1 = [NaN]*ones(sys.m,1);
        end


        diagnostics2 = optimize(Constraints2,J,options);
        if ~diagnostics2.problem
            J2 = value(J);
            u2 = value(u_hat(:,1));
        else
            J2 = Inf;
            u2 = [NaN]*ones(sys.m,1);
        end

        if diagnostics2.problem && diagnostics1.problem
           error('Infeasibility due to PE') 
        end

        if J1<J2
            u = u1;
        else
            u = u2;
        end
    else
        % PE is already satisfied. No additional constraint
        diagnostics = optimize(Constraints,J,options);
        if diagnostics.problem
           error('Infeasibility even without PE') 
        end
        u = value(u_hat(:,1));
    end
catch
   warn_flag = 1;
   warning('on');
   warning('PE tried and failed. Implementing standard RAMPC')
    diagnostics = optimize(Constraints,J,options);
    if diagnostics.problem
       error('Infeasibility even without PE') 
    end
    u = value(u_hat(:,1));
end
J_OL = value(J);
end

%% PE condition for SISO system
% 

%% Original PE condition (nonconvex LMI)
% omega_bar_PE = -cont.rho_PE*eye(sys.m*sys.n-sys.m);
% for k =0:cont.P_PE
%     phi_vec = vec(U_past(:,k+1:k+sys.n-1));
%     omega_bar_PE = omega_bar_PE + phi_vec*phi_vec';
% end   
% phi_vec0 = vec(U_past(:,1:sys.n-1));
% 
% alpha_PE = 1-phi_vec0'*(omega_bar_PE\phi_vec0);
% beta_PE = zeros(sys.m,1);
% u_phi_vec = zeros(sys.m,1);
% for j = 1:cont.P_PE
%     phi_vec = vec(U_past(:,j+1:j+sys.n-1));
%     beta_PE = beta_PE - U_past(:,j)*phi_vec'*(omega_bar_PE\phi_vec0);
%     u_phi_vec = u_phi_vec + U_past(:,j)*phi_vec';
% end
% gamma_PE = phi_vec0'*phi_vec0 - cont.rho_PE*eye(sys.m)...
%             - u_phi_vec'*(omega_bar_PE\u_phi_vec);
% u_k_var = [cont.K*xk + v_lk(:,1)];
% Constraints = [Constraints; 
%            alpha_PE*u_k_var*u_k_var' + u_k_var*beta_PE' + ...
%            beta_PE*u_k_var' + gamma_PE>=0];