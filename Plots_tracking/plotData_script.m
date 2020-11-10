%% Script to perform Robust Adaptive MPC

clc;
clear all;
yalmip clear;

addpath('..\Functions','..\')
sys = system_desc();
cont = cont_params(sys);
%% Define controller parameters

cont.theta_hat = [0.5;0.5];
cont.x_hat_k = sys.x0;
cont.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont.theta_hat,[1,1,sys.p])),3);
cont.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont.theta_hat,[1,1,sys.p])),3);

% Exploration: number of predictions
cont.nPred_theta = 1;
%% Define simulation parameters


% reference trajectory
x_s = [0.0  2.5  
       0.0  -2.0  ];
n_ref = size(x_s,2);

T_ref = [7; 18;];

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
%% Setup simulations

PAMPC_sims = 1;
DAMPC_sims = 1:7;

PE_rhos = [1,3,5,7, 8, 10, 20, 30, 40, 50, 70 80,  100, 200, 300, 400, 500,800]*1e-2;

PE_sims = 1:length(PE_rhos);

ind = PAMPC_sims + DAMPC_sims(end) + PE_sims(end);

for sims = 1:ind
    % initialize states and inputs

    clear functions;
    clear global;
    
    x = NaN*ones(sys.n,Tsim+1);
    u = NaN*ones(sys.m,Tsim);
    u_std = NaN*ones(sys.m,Tsim);

    % reset parameters
    cont.h_theta_k = cont.h_theta_0;
    cont.theta_hat = [0.5;0.5];
    
    J = 0;
    % regressors
    Dk = zeros(sys.n*cont.blk,sys.p);
    dk = zeros(sys.n*cont.blk,1);  

    % initial condition
    x(:,1) = sys.x0;

    % model initialization
    true_sys = model(sys,x(:,1));

    presolve = 0;
    if presolve
        optProb1 = controller_pre(sys,cont);
    %     optProb2 = controller_expl_pre(sys,cont);
    end

    x_ref = horzcat(ref(1:Tsim+1).x_s);

    if sims>PAMPC_sims && sims<=(PAMPC_sims+DAMPC_sims(end))
        cont.nPred_X = sims-1;
    elseif sims>PAMPC_sims+DAMPC_sims
        PE = 2;
        cont.P_PE = sys.n + 1;
        u_past_sim = [0 0 0 0
                      0 0 0 0];
        cont.rho_PE = PE_rhos(sims-PAMPC_sims-DAMPC_sims(end));
    end
    
%% Simulate system
    for k = 1:Tsim
        % update feasible set and parameter estimate
        cont = updateParameters(sys,cont,x(:,k),Dk,dk);
        h_theta_k_mat(:,k) = cont.h_theta_k;
        
        if k == 1
        rng(14,'twister');
        end

        % calculate control input   
        if sims>PAMPC_sims && sims<=PAMPC_sims+DAMPC_sims(end)
           [u_P(:,k),J_OL_P] = controller(sys,cont,x(:,k),ref(k:k+cont.N));  
        tic
           [u(:,k),J_OL(k,1)] = controller_expl(sys,cont,x(:,k),ref(k:k+cont.N)); 
           warn_flag(k,1) = 0;
           if J_OL_P<J_OL(k,1)
               warn_flag(k,1) = 1;
               warning('The optimization did gave worse solution than PAMPC');
           end
        elseif sims>PAMPC_sims+DAMPC_sims(end)
            % Build past input vector  (for PE)  
            U_past = [u(:,k-1:-1:max(1,k-sys.n-cont.P_PE)),u_past_sim(:,1:max(0,sys.n+cont.P_PE-k))];
        tic
            [u(:,k),J_OL(k,1),warn_flag(k,1)] = controller_PE(sys,cont,x(:,k),U_past,PE,ref(k:k+cont.N));
        else 
        tic
            [u(:,k),J_OL(k,1)] = controller(sys,cont,x(:,k),ref(k:k+cont.N));  
            warn_flag(k,1) = 0;
        end
        runtimes(k) = toc;
        
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
%% Save data
data(sims).x = x;
data(sims).u = u;
data(sims).J = J;
data(sims).sys = sys;
data(sims).cont = cont;
data(sims).h_theta_k_mat = h_theta_k_mat;
data(sims).x_ref = x_ref;
data(sims).runtimes = runtimes;
data(sims).warn_flag = warn_flag;
end
%%
indices.PAMPC = PAMPC_sims;
indices.DAMPC = indices.PAMPC(end)+DAMPC_sims;
indices.PE = indices.DAMPC(end)+PE_sims;
save('simulationData','data','indices');