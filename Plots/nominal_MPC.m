% script to run nominal MPC for all the systems in the parameter space

clc;
clear;

addpath('..\Functions','..\')
sys = system_desc();
cont = cont_params(sys);

Tsim = 10;
%%

% range of thetas to be considered
pts = 51;
theta1 = linspace(-0.4,0.7,pts);
theta2 = linspace(-0.5,1,pts);
J_all = NaN*ones(length(theta1));
u1_all = NaN*ones(length(theta1));
u2_all = NaN*ones(length(theta1));


tic
optProb = controller_nom(sys,cont);
toc

for t1 = 1:length(theta1)
for t2 = 1:length(theta2)
%% Simulate for each set of thetas
% initialize states and inputs
x = NaN*ones(sys.n,Tsim+1);
u = NaN*ones(sys.m,Tsim);

J = 0;

% initial condition
x(:,1) = sys.x0;

% model initialization
true_sys = model_nom(sys,[theta1(t1);theta2(t2)],x(:,1));

rng(10,'twister');
% simulate
for k = 1:Tsim
    
    % calculate control input   
    [u(:,k)] = optProb([x(:,k);theta1(t1);theta2(t2)]);   
%     runtimes(:,k) = a5.solvertime;
    
    % Apply to true system
    true_sys = true_sys.simulate(u(:,k));
    x(:,k+1) = true_sys.x;
    J = J + norm(cont.Q_L*x(:,k+1),'inf') + norm(cont.R_L*u(:,k),'inf');
    
end
%% Store relevant data

J_all(t1,t2) = J;
u1_all(t1,t2) = u(1,1);
u2_all(t1,t2) = u(2,1);

clear true_sys
end
toc
end

%% Save data
save('nom_data.mat','J_all','u1_all','u2_all')