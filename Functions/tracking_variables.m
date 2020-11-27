function [K_val,pwc_var] = tracking_variables(sys,cont,ref)
% Compute inputs and terminal sets for all values of pwc reference in ref
% The terminal controller has the form: u = K(x-x_s)+ u_s
% The terminal set has the form: X_T = H_x(x-x_s) <= alpha_T*1.
% The values of u_s and alpha_T corresponding to x_s are computed here.

Ac = NaN*ones(sys.n,sys.n,sys.p);
Bc = NaN*ones(sys.n,sys.m,sys.p);
for k  = 1:size(sys.H_theta_v,2)
    Ac(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3);
    Bc(:,:,k) = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3);       
end 


%% Compute alpha_T and u_s for each of the references

pwc_var.x_s = unique(ref','rows','stable')';
n_ref = size(pwc_var.x_s,2);

A_hat = sys.A0+sum(bsxfun(@times,sys.Ap,cont.theta_hat),3);
B_hat = sys.B0+sum(bsxfun(@times,sys.Bp,cont.theta_hat),3);

for i = 1:n_ref
    % compute input, for system corresponding to theta=0
%     pwc_var.u_s(:,i) = inv(sys.B0)*(eye(sys.n)-sys.A0)*pwc_var.x_s(:,i);

    % compute input, for system corresponding to theta=theta_hat
    pwc_var.u_s(:,i) = inv(B_hat)*(eye(sys.n)-A_hat)*pwc_var.x_s(:,i);
    
    % check constraints for the obtained control input
    chk = sys.F*pwc_var.x_s(:,i) + sys.G*pwc_var.u_s(:,i)-sys.vec_1_cons;
    if any(chk>0)
        error('Steady state input does not satisfy constraints')
    end
end
    
alpha_inv = sdpvar(n_ref,1,'full');
K_var = sdpvar(sys.m,sys.n,'full');
Cons_a = [];
for i = 1:n_ref
    for k = 1:size(sys.H_theta_v,2)
        for j = 1:size(cont.x_v,2)
           Cons_a = [Cons_a; ...
               cont.H_x*(Ac(:,:,k)+Bc(:,:,k)*K_var)*cont.x_v(:,j)+  cont.H_x*D_mult(sys,pwc_var.x_s(:,i),pwc_var.u_s(:,i))*sys.H_theta_v(:,k)*alpha_inv(i)  + cont.w_bar*alpha_inv(i)  <= cont.vec_1_x];
        end
    end

    for j = 1:size(cont.x_v,2)
       Cons_a = [Cons_a; (sys.F+sys.G*K_var)*cont.x_v(:,j)<=sys.vec_1_cons*alpha_inv(i)];
    end
end
options = sdpsettings('verbose',0);
diagnostics = optimize(Cons_a,alpha_inv,options);
if diagnostics.problem
   error(diagnostics.info) 
end

% store the value of alpha_T
pwc_var.alpha_T = value(1./alpha_inv); 
K_val = value(K_var);

end

