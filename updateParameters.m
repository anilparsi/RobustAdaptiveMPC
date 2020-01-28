function cont = updateParameters(sys,cont,xk,Dk,dk)

H_til = cont.H_theta;
h_til = cont.h_theta_k;

%% update parameter bounds
for i = 1:cont.blk
    idx1 = sys.n*(i-1)+1;
    idx2 = sys.n*i;
    
    H_til = [H_til
             -sys.H_w*Dk(idx1:idx2,:);];
         
    h_til = [h_til;   
             ones(length(sys.h_w),1) + sys.H_w*dk(idx1:idx2)];                
end
            
r = length(cont.h_theta_k);
r2 = length(h_til);

cvx_begin quiet
    variable h(r,1)
    variable Lam(r,r2)
    
    minimize(sum(h))
    
    subject to
        Lam(:)>=0;
        Lam*H_til==cont.H_theta;
        Lam*h_til<=h;
cvx_end

if ~all(cvx_status=='Solved')
   warning(cvx_status) 
   flag = 1;
   u = [];
   return;
end

% update bounds of polytope
cont.h_theta_k = h;
    
%% Estimate of theta
theta_til = cont.theta_hat + cont.mu * Dk(end-sys.n+1:end,:)'*(xk-cont.x_hat_k); 

if any(cont.H_theta*theta_til>=cont.h_theta_0)
    % theta_til outside the initial bounds, use projection
    cvx_begin 
        variable theta_hat_p(sys.p,1)
        minimize sum((theta_hat_p-theta_til).^2)
        
        subject to
            cont.H_theta*theta_hat_p <= cont.h_theta_0;
        
    cvx_end
    
    if ~all(cvx_status=='Solved')
       warning(cvx_status) 
       flag = 1;
       u = [];
       return;
    else 
       cont.theta_hat = theta_hat_p;
    end
    
else
    cont.theta_hat = theta_til;   
end

cont.A_est = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(cont.theta_hat,[1,1,3])),3);
cont.B_est = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(cont.theta_hat,[1,1,3])),3);
end