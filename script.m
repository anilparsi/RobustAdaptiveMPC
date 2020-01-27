%% Script to perform Robust Adaptive MPC

clc;
clear;

sys = system_desc();
%% Define controller parameters

% prediction horizon
cont.N = 9;

% cost matrices
cont.Q = eye(sys.n);
cont.R = eye(sys.m);
cont.P = [1.467,0.207;
          0.207,1.731];

% block size for updating h_theta_k
cont.blk = 10;

% prestabilizing gain
cont.K = [0.017 -0.41];

% parameter bounds
cont.H_theta = sys.H_theta;
cont.h_theta_k = sys.h_theta;
cont.h_theta_0 = sys.h_theta;

% Estimation parameters
cont.mu = 2; % compute properly later
cont.theta_hat = [0.5;0.5;0.5];
cont.x_hat_k = sys.x0;

% Define state tube shape : H_x*x <= vec_1_x, vertices: x_v(:,i)
cont.H_x = [eye(sys.n); -eye(sys.n)];
cont.vec_1_x = ones(length(cont.H_x),1);
cont.nHx = size(cont.H_x,1);  % u in Lorenzen(2019)
cont.x_v = [1 1; 1 -1; -1 1; -1 -1 ]';
cont.nx_v = length(cont.x_v); % number of vertices

% Calculate constants fBar and wBar
cont.f_bar = max((sys.F+sys.G*cont.K)*cont.x_v,[],2);
cont.w_bar = [0.1; 0.1; 0.1; 0.1];  % manually calculated. need to verify if disturbance changed

% Define terminal constraint: z_N|k = 0; h_T*alpha_N|k <= 1;
cont.h_T = 1;

% Set cvx_solver
cvx_solver gurobi
%% Define simulation parameters

Tsim = 10;

% initialize states and inputs
x = NaN*ones(sys.n,Tsim);
u = NaN*ones(sys.m,Tsim);

theta_hats = NaN*ones(sys.p,Tsim);


% regressors
Dk = zeros(sys.n*cont.blk,sys.p);
dk = zeros(sys.n*cont.blk,1);  

% initial condition
x(:,1) = sys.x0;

% model initialization
true_sys = model(sys,x(:,1));

% simulate
for k = 1:Tsim
    % update feasible set and parameter estimate
    cont = updateParameters(sys,cont,x(:,k),Dk,dk);
    theta_hats(:,k) = cont.theta_hat;
    
    % calculate control input
    u(:,k)= controller(sys,cont,x(:,k));
    
    % update state estimate
    cont.x_hat_k = cont.A_est*x(:,k)+cont.B_est*u(:,k);
    
    % Apply to true system
    true_sys = true_sys.simulate(u(:,k));
    x(:,k+1) = true_sys.x;

    
    % update regressors
        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x(:,k) + sys.Bp(:,:,i)*u(:,k);
        end    
        new_dk = -x(:,k+1)+sys.A0*x(:,k)+sys.B0*u(:,k);

        Dk = [Dk(sys.n+1:end,:);new_Dk];
        dk = [dk(sys.n+1:end);new_dk];
        
end

%%
plot(x(1,:),x(2,:))