%% Script to perform Robust Adaptive MPC

clc;
clear;

sys = system_desc();
%% Define controller parameters

% prediction horizon
cont.N = 9;


% cost matrices
cont.Q = eye(2);
cont.R = eye(1);
cont.P = [1.467,0.207;
    0.207,1.731];

% block size for updating h_theta_k
cont.blk = 10;

% prestabilizing gain
cont.K = [0.017 -0.41];

% parameter bounds
cont.H_theta = sys.H_theta;
cont.h_theta_k = sys.h_theta;

cont = boundedComplexity(cont);

%% Define simulation parameters

Tsim = 100;

% initialize states and inputs
x = NaN*ones(sys.n,Tsim);
u = NaN*ones(sys.m,Tsim);


% regressors
Dk = zeros(sys.n*cont.blk,sys.p);
dk = zeros(sys.n*cont.blk,1);  

% initial condition
x(:,1) = [2;3];

% model initialization
true_sys = model(sys,x(:,1));

% cont.H_theta(7:end,:) = [];
% cont.h_theta_k(7:end,:) = [];

for t = 1:Tsim
    [u(:,t),cont] = controller(sys,cont,x(:,t),Dk,dk);
    
    true_sys = true_sys.simulate(u(:,t));
    x(:,t+1) = true_sys.x;

    % update regressors
        Dk(1:sys.n,:) = [];
        dk(1:sys.n) = [];

        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x(:,t) + sys.Bp(:,:,i)*u(:,t);
        end    
        new_dk = -x(:,t+1)+sys.A0*x(:,t)+sys.B0*u(:,t);

        Dk = [Dk;new_Dk];
        dk = [dk;new_dk];
end