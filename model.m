classdef model
properties (Constant)
    theta_true = [0.8; 0.2; -0.5];
end
properties
    A_true 
    B_true 
    w_bound
    x
end


methods
    % Constructor
    function obj = model(sys,x0)
        obj.A_true = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(obj.theta_true,[1,1,3])),3);
        obj.B_true = sys.B0+ sum(bsxfun(@times,sys.Bp,reshape(obj.theta_true,[1,1,3])),3);
        
        obj.w_bound = sys.w_bound;
        obj.x = x0;
    end
    
    % simulate
    function obj = simulate(obj,u)
        % works only for uniform disturbance/noise bounds
        w = -obj.w_bound + 2*obj.w_bound*rand(length(obj.x),1);        
        
        obj.x = obj.A_true*obj.x + obj.B_true*u+w;
    end
    
end

end