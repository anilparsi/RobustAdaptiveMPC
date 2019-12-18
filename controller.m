function [u,cont]= controller(sys,cont,xt,Dk,dk)

H_til = cont.H_theta;
h_til = cont.h_theta_k;

% update parameter bounds
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

cvx_begin %quiet
cvx_solver gurobi
    variable h(r,1)
    variable Lam(r,r2)
    
    minimize(sum(h))
    
    subject to
        Lam(:)>=0;
        Lam*H_til==cont.H_theta;
        Lam*h_til<=h;
cvx_end

if (isnan(sum(h)))
   warning(cvx_status) 
   flag = 1;
   u = [];
   return;
end

% update bounds of polytope
cont.h_theta_k = h;
    
u = rand();
end

