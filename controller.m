function u = controller(sys,cont,xk)

%% setup optimization problem

% Declare independent variables
% l goes from 1:N
% j goes from 1:nx_v
z_lk = sdpvar(sys.n,cont.N+1,'full');
v_lk = sdpvar(sys.m,cont.N,'full');
alpha_lk = sdpvar(1,cont.N+1,'full');
Lambda_jlk = sdpvar(cont.nHx,sys.nHtheta,cont.nx_v,cont.N,'full');

% define cost variables
x_hat = xk;
u_hat = [];
J = 0;
for l =1:cont.N
    u_hat = [u_hat  cont.K*x_hat(:,l) + v_lk(:,l)] ;
    x_hat = [x_hat  cont.A_est*x_hat(:,l) + cont.B_est*u_hat(:,l)];
    J = J +  x_hat(:,l)'*cont.Q*x_hat(:,l) + u_hat(:,l)'*cont.R*u_hat(:,l);
end
J = J + x_hat(:,cont.N+1)'*cont.P*x_hat(:,cont.N+1);

% define state dependent variables
x_jlk = [];
u_jlk = [];
d_jlk = [];
D_jlk = [];
for l = 1:cont.N    
    x_jlk = cat(3,x_jlk,repmat(z_lk(:,l),1,cont.nx_v)+alpha_lk(:,l)*cont.x_v);
    u_jlk = cat(3,u_jlk,repmat(v_lk(:,l),1,cont.nx_v)+cont.K*x_jlk(:,:,l));
    
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
         (sys.F+sys.G*cont.K)*z_lk(:,l) + sys.G*v_lk(:,l) + alpha_lk(:,l)*cont.f_bar <= ones(sys.nc,1)];         
     
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
            
options = sdpsettings('solver','gurobi','verbose',0,'debug',0);

diagnostics = optimize(Constraints,J,options);
if diagnostics.problem
   error(diagnostics.info);
end
u = value(u_hat(:,1));


end

