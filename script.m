%% Script to perform Robust Adaptive MPC

clc;
clear;

sys = system_desc();
%% Define controller parameters

% prediction horizon
cont.N = 8;

% cost matrices
cont.Q = eye(sys.n);
cont.R = eye(sys.m);
cont.P = [2.8  0.56;
          0.56 2.5];
cont.Q_L = chol(cont.Q);
cont.R_L = chol(cont.R);
      
      
% block size for updating h_theta_k
cont.blk = 10;

% Define state tube shape : H_x*x <= vec_1_x, vertices: x_v(:,i)
cont.H_x = [eye(sys.n); -eye(sys.n)];
cont.vec_1_x = ones(length(cont.H_x),1);
cont.nHx = size(cont.H_x,1);  % denoted u in Lorenzen(2019)
cont.x_v = [1 1; 1 -1; -1 1; -1 -1 ]';
cont.nx_v = length(cont.x_v); % number of vertices

% prestabilizing gain
addpath('Functions')
cont.w_bar = compute_wbar(sys,cont);
[cont.K,cont.alpha_bar,cont.alpha_min] = prestab_controller(sys,cont);

% parameter bounds
cont.H_theta = sys.H_theta;
cont.h_theta_k = sys.h_theta;
cont.h_theta_0 = sys.h_theta;

% Estimation parameters
cont.mu = 2; 
cont.theta_hat = [0.5;0.5;];
cont.x_hat_k = sys.x0;
cont.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont.theta_hat,[1,1,sys.p])),3);
cont.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont.theta_hat,[1,1,sys.p])),3);


% Define terminal constraint: z_N|k = 0; h_T*alpha_N|k <= 1;
cont.f_bar = max((sys.F+sys.G*cont.K)*cont.x_v,[],2);

% Exploration: number of predictions
cont.nPred_theta = 1;
cont.nPred_X = 7;
%% Define simulation parameters

Tsim = 10;

% initialize states and inputs
x = NaN*ones(sys.n,Tsim+1);
u = NaN*ones(sys.m,Tsim);
u_std = NaN*ones(sys.m,Tsim);

theta_hats = NaN*ones(sys.p,Tsim);
J = 0;
% regressors
Dk = zeros(sys.n*cont.blk,sys.p);
dk = zeros(sys.n*cont.blk,1);  

% initial condition
x(:,1) = sys.x0;

% model initialization
true_sys = model(sys,x(:,1));

tic
presolve = 1;
if presolve
%     optProb1 = controller_pre(sys,cont);
    optProb2 = controller_expl_pre(sys,cont);
end
toc

rng(10,'twister');
% simulate
for k = 1:Tsim
    if any(k == [15, 30, 50])
        figure; plotregion(-cont.H_theta,-cont.h_theta_k);
        xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    end
    % update feasible set and parameter estimate
    cont = updateParameters(sys,cont,x(:,k),Dk,dk);
    theta_hats(:,k) = cont.theta_hat;
    
    % calculate control input   
    tic
    if presolve
%         [u(:,k),a1,~,~,~,a5] = optProb1([x(:,k);cont.h_theta_k]);
        [u(:,k),a1,~,~,~,a5] = optProb2([x(:,k);cont.h_theta_k;cont.theta_hat]);
        if a1
            error(a5.infostr)
        end
    else
        u_std(:,k) = controller(sys,cont,x(:,k));   
        u(:,k) = controller_expl(sys,cont,x(:,k));    
    end
    toc
    
    % update state estimate
    cont.x_hat_k = cont.A_est*x(:,k)+cont.B_est*u(:,k);
    
    % Apply to true system
    true_sys = true_sys.simulate(u(:,k));
    x(:,k+1) = true_sys.x;
    J = J + norm(cont.Q_L*x(:,k+1),'inf') + norm(cont.R_L*u(:,k),'inf');
    
    % update regressors
        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x(:,k) + sys.Bp(:,:,i)*u(:,k);
        end    
        new_dk = -x(:,k+1)+sys.A0*x(:,k)+sys.B0*u(:,k);

        Dk = [Dk(sys.n+1:end,:);new_Dk];
        dk = [dk(sys.n+1:end);new_dk];
        
end
toc
%%
f = 60;
h = figure(f); f=f+1;
subplot(411)
hold on;
plot(0:Tsim,x(1,:),'-*')
ylabel('$x_1$')

subplot(412)
hold on;
plot(0:Tsim,x(2,:),'-*')
ylabel('$x_2$')

subplot(413)
hold on;
plot(0:Tsim-1,u(1,:),'-*')
xlim([0 Tsim])
ylabel('$u_1$')

subplot(414)
hold on;
plot(0:Tsim-1,u(2,:),'-*')
xlim([0 Tsim])
ylabel('$u_2$')
% h = figure(f);f=f+1;
% plotregion(-cont.H_theta,-cont.h_theta_k);
% xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);