function optProb = controller_nom(sys,cont)

%% setup optimization problem

x_pre = sdpvar(sys.n,1,'full');
theta_pre = sdpvar(sys.p,1,'full');

%% NOTE: HARDCODED FOR 2 dimensional parameters
A_true = sys.A0+ sys.Ap(:,:,1)*theta_pre(1)+ sys.Ap(:,:,2)*theta_pre(2);
B_true = sys.B0+ sys.Bp(:,:,1)*theta_pre(1)+ sys.Bp(:,:,2)*theta_pre(2);
        
% Declare independent variables
% l goes from 1:N
% j goes from 1:nx_v
x_lk = sdpvar(sys.n,cont.N+1,'full');
v_lk = sdpvar(sys.n,cont.N,'full');
Constraints = [x_lk(:,1) == x_pre];
for l = 1:cont.N
    % constraints
     Constraints = [Constraints, ...
         (sys.F+sys.G*cont.K)*x_lk(:,l) + sys.G*v_lk(:,l) <= ones(sys.nc,1)];         
    % dynamics
     Constraints = [Constraints, ...
         x_lk(:,l+1) == (A_true+ B_true *cont.K)*x_lk(:,l) + B_true *v_lk(:,l)];
end

% define terminal constraints
Constraints = [Constraints, ...
              cont.H_x *x_lk(:,cont.N+1) <= cont.vec_1_x*cont.alpha_bar;];
            
options = sdpsettings('solver','gurobi','verbose',0,'debug',0);

%% Cost function definition

J = 0;

for l = 1:cont.N
    J = J + (norm(cont.Q_L*x_lk(:,l),'inf') + norm(cont.R_L*(cont.K*x_lk(:,l)+v_lk(:,l)),'inf'));
end
J = J + (norm(cont.Q_L*x_lk(:,l),'inf') + norm(cont.R_L*(cont.K*x_lk(:,l)),'inf'));

optProb = optimizer(Constraints,J,options,[x_pre;theta_pre],cont.K*x_pre+v_lk(:,1));


end

