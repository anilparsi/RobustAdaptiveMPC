function K_val = prestab_controller(sys,cont)
% function to find a robust control gain matrix

warning('Using hardcoded vertices for H_theta. Check their validity')

A_corners = NaN*ones(sys.n,sys.n,sys.p);
B_corners = NaN*ones(sys.n,sys.m,sys.p);
for k  = 1:size(sys.H_theta_v,1)
    A_corners(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(k,:),[1,1,3])),3);
    B_corners(:,:,k) = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(sys.H_theta_v(k,:),[1,1,3])),3);       
end


K = sdpvar(sys.m,sys.n,'full');
Constraints = [];
    for k = 1:size(sys.H_theta_v,1)
     Constraints = [Constraints, ...        
    [cont.P-cont.Q          K'          (A_corners(:,:,k)+B_corners(:,:,k)*K)'
            K                         cont.R^-1             zeros(sys.m,sys.n)
A_corners(:,:,k)+B_corners(:,:,k)*K   zeros(sys.n,sys.m)    inv(cont.P)
        ] >= 0];
    end
    
options = sdpsettings('verbose',1);
diagnostics = optimize(Constraints,[],options);
K_val = value(K);

for k = 1:size(sys.H_theta_v,1)
    eigs = eig(A_corners(:,:,k)+B_corners(:,:,k)*K_val);
    if any(eigs>1)
        error('Gain computed is wrong')     
    end
end


end