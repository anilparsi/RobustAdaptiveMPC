%% Script to generate plots for RAMPC

clc;
clear;

addpath('..\Functions','..\')
sys = system_desc();
cont_P = cont_params(sys);

Tsim = 10;

%% Simulate PAMPC

% initialize states and inputs
x_P = NaN*ones(sys.n,Tsim+1);
u_P = NaN*ones(sys.m,Tsim);
runtimes_P = NaN*ones(Tsim,1);

J_P = 0;
% regressors
Dk_P = zeros(sys.n*cont_P.blk,sys.p);
dk_P = zeros(sys.n*cont_P.blk,1);  


% initial condition
x_P(:,1) = sys.x0;

% model initialization
true_sys = model(sys,x_P(:,1));

% Estimation parameters % Unused in PAMPC (needed in sim only)
cont_P.mu = 2; 
cont_P.theta_hat = [0.5;0.5;];
cont_P.x_hat_k = sys.x0;
cont_P.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont_P.theta_hat,[1,1,sys.p])),3);
cont_P.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont_P.theta_hat,[1,1,sys.p])),3);

optProb_P = controller_pre(sys,cont_P);

rng(10,'twister');
% simulate
for k = 1:Tsim
    % update feasible set and parameter estimate
    cont_P = updateParameters(sys,cont_P,x_P(:,k),Dk_P,dk_P);
    
    % calculate control input   
    [u_P(:,k),a1,~,~,~,a5] = optProb_P([x_P(:,k);cont_P.h_theta_k]);   
    runtimes_P(:,k) = a5.solvertime;
    
    % Apply to true system
    true_sys = true_sys.simulate(u_P(:,k));
    x_P(:,k+1) = true_sys.x;
    J_P = J_P + norm(cont_P.Q_L*x_P(:,k+1),'inf') + norm(cont_P.R_L*u_P(:,k),'inf');
    
    % update regressors
        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x_P(:,k) + sys.Bp(:,:,i)*u_P(:,k);
        end    
        new_dk = -x_P(:,k+1)+sys.A0*x_P(:,k)+sys.B0*u_P(:,k);

        Dk_P = [Dk_P(sys.n+1:end,:);new_Dk];
        dk_P = [dk_P(sys.n+1:end);new_dk];
        
end

%% Simulate DAMPC 1
cont_D1 = cont_params(sys);

% initialize states and inputs
x_D1 = NaN*ones(sys.n,Tsim+1);
u_D1 = NaN*ones(sys.m,Tsim);
runtimes_D1 = NaN*ones(Tsim,1);

J_D1 = 0;
% regressors
Dk_D1 = zeros(sys.n*cont_D1.blk,sys.p);
dk_D1 = zeros(sys.n*cont_D1.blk,1);  

% initial condition
x_D1(:,1) = sys.x0;

% model initialization
true_sys = model(sys,x_D1(:,1));

% Estimation parameters
cont_D1.mu = 2; 
cont_D1.theta_hat = [0.5;0.5;];
cont_D1.x_hat_k = sys.x0;
cont_D1.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont_D1.theta_hat,[1,1,sys.p])),3);
cont_D1.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont_D1.theta_hat,[1,1,sys.p])),3);

% Exploration: number of predictions
cont_D1.nPred_theta = 1;
cont_D1.nPred_X = 2;

optProb_D1 = controller_expl_pre(sys,cont_D1);

rng(10,'twister');
% simulate
for k = 1:Tsim
    % update feasible set and parameter estimate
    cont_D1 = updateParameters(sys,cont_D1,x_D1(:,k),Dk_D1,dk_D1);
    
    % calculate control input   
    [u_D1(:,k),a1,~,~,~,a5] = optProb_D1([x_D1(:,k);cont_D1.h_theta_k;cont_D1.theta_hat]);   
    runtimes_D1(k) = a5.solvertime;
    
    % Apply to true system
    true_sys = true_sys.simulate(u_D1(:,k));
    x_D1(:,k+1) = true_sys.x;
    J_D1 = J_D1 + norm(cont_D1.Q_L*x_D1(:,k+1),'inf') + norm(cont_D1.R_L*u_D1(:,k),'inf');
    
    % update regressors
        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x_D1(:,k) + sys.Bp(:,:,i)*u_D1(:,k);
        end    
        new_dk = -x_D1(:,k+1)+sys.A0*x_D1(:,k)+sys.B0*u_D1(:,k);

        Dk_D1 = [Dk_D1(sys.n+1:end,:);new_Dk];
        dk_D1 = [dk_D1(sys.n+1:end);new_dk];
end

%% Simulate DAMPC 2

cont_D2 = cont_params(sys);

% initialize states and inputs
x_D2 = NaN*ones(sys.n,Tsim+1);
u_D2 = NaN*ones(sys.m,Tsim);
runtimes_D2 = NaN*ones(Tsim,1);

J_D2 = 0;
% regressors
Dk_D2 = zeros(sys.n*cont_D2.blk,sys.p);
dk_D2 = zeros(sys.n*cont_D2.blk,1);  

% initial condition
x_D2(:,1) = sys.x0;

% model initialization
true_sys = model(sys,x_D2(:,1));

% Estimation parameters
cont_D2.mu = 2; 
cont_D2.theta_hat = [0.5;0.5;];
cont_D2.x_hat_k = sys.x0;
cont_D2.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont_D2.theta_hat,[1,1,sys.p])),3);
cont_D2.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont_D2.theta_hat,[1,1,sys.p])),3);

% Exploration: number of predictions
cont_D2.nPred_theta = 1;
cont_D2.nPred_X = 5;

optProb_D2 = controller_expl_pre(sys,cont_D2);

rng(10,'twister');
% simulate
for k = 1:Tsim
    % update feasible set and parameter estimate
    cont_D2 = updateParameters(sys,cont_D2,x_D2(:,k),Dk_D2,dk_D2);
    
    % calculate control input   
    [u_D2(:,k),a1,~,~,~,a5] = optProb_D2([x_D2(:,k);cont_D2.h_theta_k;cont_D2.theta_hat]);   
    runtimes_D2(k) = a5.solvertime;
    
    % Apply to true system
    true_sys = true_sys.simulate(u_D2(:,k));
    x_D2(:,k+1) = true_sys.x;
    J_D2 = J_D2 + norm(cont_D2.Q_L*x_D2(:,k+1),'inf') + norm(cont_D2.R_L*u_D2(:,k),'inf');
    
    % update regressors
        new_Dk = zeros(sys.n,sys.p);    
        for i= 1:sys.p
            new_Dk(:,i) = sys.Ap(:,:,i)*x_D2(:,k) + sys.Bp(:,:,i)*u_D2(:,k);
        end    
        new_dk = -x_D2(:,k+1)+sys.A0*x_D2(:,k)+sys.B0*u_D2(:,k);
    
        Dk_D2 = [Dk_D2(sys.n+1:end,:);new_Dk];
        dk_D2 = [dk_D2(sys.n+1:end);new_dk];
        
end

%% PLOTS START
%% Plot the trajectories

set(0,'DefaultFigureWindowStyle','normal','defaultLineLineWidth',1.5,'DefaultLineMarkerSize',3)

f = 80; 
h = figure(f);clf; f = f+1; 
h.Units = 'centimeters';
h.WindowStyle = 'normal';
h.Position = [2 2 8 10];
ha = tight_subplot(4,1,[0.03 0.05],[.09 .01],[.1 .03]);

i = 1;
axes(ha(i));hold on;
plot(0:Tsim,x_P(1,:),'b-^')
plot(0:Tsim,x_D1(1,:),'r-.^')
plot(0:Tsim,x_D2(1,:),'k:^')
ylabel(ha(i),'$x_1$')
ylim(ha(i),[-0.2 1.5])
yticks(ha(i),0:1)
legend('PAMPC','DAMPC$_2$','DAMPC$_5$',...
      'Orientation','Vertical','Location','Northeast')

i = i+1;
axes(ha(i));hold on;
plot(0:Tsim,x_P(2,:),'b-^')
plot(0:Tsim,x_D1(2,:),'r-.^')
plot(0:Tsim,x_D2(2,:),'k:^')
ylabel(ha(i),'$x_2$')
ylim(ha(i),[-0.2 1.5])
yticks(ha(i),0:1)

i = i+1;
axes(ha(3));hold on;
plot(0:Tsim-1,-0.5*ones(Tsim,1),'m--')
plot(0:Tsim-1,u_P(1,:),'b-^')
plot(0:Tsim-1,u_D1(1,:),'r-.^')
plot(0:Tsim-1,u_D2(1,:),'k:^')
xlim([0 Tsim])
ylabel(ha(i),'$u_1$')
yticks(ha(i),-1:0)
ylim(ha(i),[-1.2 0.2])
legend('Input constraint','Location','Southeast')

i = i+1;
axes(ha(i));hold on;
plot(0:Tsim-1,u_P(2,:),'b-^')
plot(0:Tsim-1,u_D1(2,:),'r-.^')
plot(0:Tsim-1,u_D2(2,:),'k:^')
xlim([0 Tsim])
ylabel(ha(i),'$u_2$')
yticks(ha(i),-1:0)
ylim(ha(i),[-1.2 0.2])

yticklabels(ha,'auto')
xticklabels(ha(4),'auto')
xlabel(ha(4),'Timesteps')

export_fig trajCL.eps
% saveas(f-1,'trajCL','epsc')

%% Plot parameter set after 10 time steps

h = figure(f);clf; f = f+1; 
h.Units = 'centimeters';
h.WindowStyle = 'normal';
h.Position = [2 2 8 7];

hold on;
plotregionLine(-cont_P.H_theta,-cont_P.h_theta_k,[],[],'b',[],[],'-');
plotregionLine(-cont_D1.H_theta,-cont_D1.h_theta_k,[],[],'r',[],[],'-.');
plotregionLine(-cont_D2.H_theta,-cont_D2.h_theta_k,[],[],'k',[],[],'--');
plot(0.9,0.3,'m^')
xlim([-1 1]);ylim([-1 1]);
xlabel('$\theta_1$')
ylabel('$\theta_2$')
legend('PAMPC','DAMPC$_2$','DAMPC$_5$',...
      'Orientation','Vertical','Location','southeast')
  
export_fig parameterSet.eps