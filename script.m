%% Script to perform Robust Adaptive MPC

clc;
clear all;
yalmip clear;
sys = system_desc();
%% Define controller parameters

% prediction horizon
cont.N = 8;

% cost matrices
cont.Q = eye(sys.n);
cont.R = 1*eye(sys.m);
cont.P = [2.8  0.56;
          0.56 2.5];
cont.Q_L = cont.Q;
cont.R_L = cont.R;
      
      
% block size for updating h_theta_k
cont.blk = 5;

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
cont.theta_hat = [0.5;0.5];
cont.x_hat_k = sys.x0;
cont.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont.theta_hat,[1,1,sys.p])),3);
cont.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont.theta_hat,[1,1,sys.p])),3);

% Exploration: number of predictions
cont.nPred_theta = 1;
cont.nPred_X = 7;
%% Define simulation parameters


% reference trajectory
x_s = [0.0  2.5  
       0.0  -2.0  ]*0.5;
n_ref = size(x_s,2);

T_ref = [5; 8];

x_ref = [];
for i = 1:n_ref
    x_ref = [x_ref repmat(x_s(:,i),1,T_ref(i))];
end
Tsim = size(x_ref,2);

% find steady state inputs and terminal sets
[cont.K,pwc_var] = tracking_variables(sys,cont,x_ref);
cont.f_bar = max((sys.F+sys.G*cont.K)*cont.x_v,[],2);
u_ref = [];
alpha_T_ref = [];
for i = 1:n_ref
    u_ref = [u_ref repmat(pwc_var.u_s(:,i),1,T_ref(i))];
    alpha_T_ref = [alpha_T_ref repmat(pwc_var.alpha_T(i),1,T_ref(i))];
end

% stack reference variables
for k = 1:Tsim
    ref(k).x_s = x_ref(:,k);
    ref(k).u_s = u_ref(:,k);
    ref(k).alpha_T = alpha_T_ref(k);
end
for i = 1:cont.N
    ref(Tsim+i) = ref(end);
end

% create unique variables for reference
% % x_s already exists
% u_s = unique(horzcat(ref.u_s)','rows','stable')';
% alpha_s = unique(horzcat(ref.alpha_T),'rows','stable');

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
presolve = 0;
if presolve
    optProb1 = controller_pre(sys,cont);
%     optProb2 = controller_expl_pre(sys,cont);
end
toc

PE = 2;
cont.rho_PE = 0.8;
cont.P_PE = sys.n + 1;
u_past_sim = [0 0 0 0
              0 0 0 0];


x_ref = horzcat(ref(1:Tsim+1).x_s);

% simulate
for k = 1:Tsim
    % update feasible set and parameter estimate
    cont = updateParameters(sys,cont,x(:,k),Dk,dk);
%     theta_hats(:,k) = cont.theta_hat;

    
% 
%%
    
    % Build past input vector  (for PE)  
    U_past = [u(:,k-1:-1:max(1,k-sys.n-cont.P_PE)),u_past_sim(:,1:max(0,sys.n+cont.P_PE-k))];
   
    if k == 1
    rng(24,'twister');
    end
    
    % calculate control input   
    tic
    if presolve
%         [u(:,k),a1,~,~,~,a5] = optProb1([x(:,k);cont.h_theta_k;vertcat(ref(k:k+cont.N).x_s);ref(k+cont.N).u_s;ref(k+cont.N).alpha_T]);
        [u(:,k),a1,~,~,~,a5] = optProb2([x(:,k);cont.h_theta_k;cont.theta_hat]);
        if a1
            error(a5.infostr)
        end
    else
%         [u(:,k),J_OL(1,k)] = controller(sys,cont,x(:,k),ref(k:k+cont.N));  
%         [u(:,k),J_OL(1,k),warn_flag] = controller_PE(sys,cont,x(:,k),U_past,PE,ref(k:k+cont.N));
        [u(:,k),J_OL(1,k)] = controller_expl(sys,cont,x(:,k),ref(k:k+cont.N));    
    end
    toc
    
    % update state estimate
    cont.x_hat_k = cont.A_est*x(:,k)+cont.B_est*u(:,k);
    
    % Apply to true system
    true_sys = true_sys.simulate(u(:,k));
    x(:,k+1) = true_sys.x;
    J = J + norm(cont.Q_L*(x(:,k+1)-ref(k+1).x_s),'inf') + norm(cont.R_L*u(:,k),'inf');
    
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

%% Plot 

f = 70;
h = figure(f);  f=f+1;
%clf;
subplot(4,2,1)
hold on;
plot(0:Tsim,x_ref(1,:),'k--')
plot(0:Tsim,x(1,:),'-*')
xlim([0 Tsim])
ylim([-1 3])
ylabel('$x_1$')

subplot(4,2,3)
hold on;
plot(0:Tsim,x_ref(2,:),'k--')
plot(0:Tsim,x(2,:),'-*')
xlim([0 Tsim])
ylim([-3 1])
ylabel('$x_2$')

subplot(4,2,5)
hold on;
plot(0:Tsim-1,u(1,:),'-*')
xlim([0 Tsim])
ylim([-4 4])
ylabel('$u_1$')

subplot(4,2,7)
hold on;
plot(0:Tsim-1,u(2,:),'-*')
xlim([0 Tsim])
ylim([-4 4])
ylabel('$u_2$')

subplot(4,2,[2,4,6,8]); 
hold on; 
plotregion(-cont.H_theta,-cont.h_theta_k);
xlabel('$\theta_1$')
ylabel('$\theta_2$')
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

