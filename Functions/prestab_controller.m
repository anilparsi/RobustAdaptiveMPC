function K_val = prestab_controller(sys,cont,check)
% function to find a robust control gain matrix
% Using method from Kothare (1996)

warning('Using hardcoded vertices for H_theta. Check their validity')

Ac = NaN*ones(sys.n,sys.n,sys.p);
Bc = NaN*ones(sys.n,sys.m,sys.p);
for k  = 1:size(sys.H_theta_v,1)
    Ac(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(k,:),[1,1,3])),3);
    Bc(:,:,k) = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(sys.H_theta_v(k,:),[1,1,3])),3);       
end

if ~isempty(check) && check
    % Just check the current gain
    for k = 1:size(sys.H_theta_v,1)
        eigs = abs(eig(Ac(:,:,k)+Bc(:,:,k)*cont.K));
        if any(eigs>1)
            error('Current gain is invalid')                 
        end        
    end
    fprintf('Current gain is valid')    
    K_val = [];
else
    % find a new control gain
    Y = sdpvar(sys.m,sys.n,'full');
    Pinv = sdpvar(sys.n,sys.n); % symmetric
    gamma = sdpvar();
    Constraints = [Pinv>=0, ...
                   [1 sys.x0'
                   sys.x0 Pinv]>=0
                  ];
        for k = 1:size(sys.H_theta_v,1)
         Constraints = [Constraints, ...        
            [Pinv   Pinv*Ac(:,:,k)'+Y'*Bc(:,:,k)'   Pinv*cont.Q_L   Y'*cont.R_L
             Ac(:,:,k)*Pinv+Bc(:,:,k)*Y     Pinv    zeros(sys.n)    zeros(sys.n,sys.m)
             cont.Q_L*Pinv      zeros(sys.n)        gamma*eye(sys.n) zeros(sys.n,sys.m)
             cont.R_L*Y      zeros(sys.m,sys.n)     zeros(sys.m,sys.n) gamma*eye(sys.m)
            ] >= 0];
        end

    options = sdpsettings('verbose',0);
    diagnostics = optimize(Constraints,gamma,options);
    if diagnostics.problem
       error(diagnostics.info) 
    end
    K_val = value(Y)/value(Pinv);

    for k = 1:size(sys.H_theta_v,1)
        eigs = abs(eig(Ac(:,:,k)+Bc(:,:,k)*K_val));
        if any(eigs>1)
            error('Gain computed is wrong')     
        end
    end
end

end