function [K_val,alpha_bar,alpha_min] = prestab_controller(sys,cont)
% function to find a robust control gain matrix
% Computes K such that it is prestabilizing, and that the terminal set is
% RPI and of the form alpha_bar*X0. Any set smaller than this will also
% work as a terminal set. 

warning('Using hardcoded vertices for H_theta. Check their validity')

Ac = NaN*ones(sys.n,sys.n,sys.p);
Bc = NaN*ones(sys.n,sys.m,sys.p);
for k  = 1:size(sys.H_theta_v,2)
    Ac(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3);
    Bc(:,:,k) = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3);       
end 


%% Finding a control gain to ensure RPI of X0
K_var = sdpvar(sys.m,sys.n);
alpha_inv = sdpvar(1);
Cons_K = [];

for k = 1:size(sys.H_theta_v,2)
    for j = 1:size(cont.x_v,2)
       Cons_K = [Cons_K; max(cont.H_x*(Ac(:,:,k)+Bc(:,:,k)*K_var)*cont.x_v(:,j))<= cont.vec_1_x - cont.w_bar*alpha_inv; ];
    end
end
for j = 1:size(cont.x_v,2)
   Cons_K = [Cons_K; (sys.F+sys.G*K_var)*cont.x_v(:,j)<=sys.vec_1_cons*alpha_inv];
end

options = sdpsettings('verbose',0);
diagnostics = optimize(Cons_K,alpha_inv,options);
if diagnostics.problem
   error(diagnostics.info) 
end
%%
% Ensure that the gain is stabilizing for all plants
K_val = value(K_var);
alpha_bar = value(1/alpha_inv);

for k = 1:size(sys.H_theta_v,2)
    eigs = abs(eig(Ac(:,:,k)+Bc(:,:,k)*K_val));
    if any(eigs>1)
        error('Gain computed is not stabilizing')     
    end
end

% Compute minimum value of alpha
alpha_min = value(alpha_bar);
for k = 1:size(sys.H_theta_v,2)
    for j = 1:size(cont.x_v,2)
        alpha_min = min(alpha_min,min(cont.w_bar./(cont.vec_1_x-cont.H_x*(Ac(:,:,k)+Bc(:,:,k)*K_val)*cont.x_v(:,j))));
    end
end

end

%%  LMI with P and K included together
% [Pinv   Pinv*Ac(:,:,k)'+Y'*Bc(:,:,k)'   Pinv*cont.Q_L   Y'*cont.R_L
%  Ac(:,:,k)*Pinv+Bc(:,:,k)*Y     Pinv    zeros(sys.n)    zeros(sys.n,sys.m)
%  cont.Q_L*Pinv      zeros(sys.n)        gamma*eye(sys.n) zeros(sys.n,sys.m)
%  cont.R_L*Y      zeros(sys.m,sys.n)     zeros(sys.m,sys.n) gamma*eye(sys.m)
% ] >= 0];