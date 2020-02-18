function X0 = state_tube(sys,cont)
% Construct initial tube shape for the system

% We start with a symmetric shape: ||x||<=1, and then construct output
% admissible sets for 2-3 time steps. This should improve performance over
% the initial set, without a complicated set description

warning('Function hardcoded for n = 2');
assert(sys.n==2);

% initial constraint set
A_nom = sys.A0;%+sys.B0*cont.K;

nSteps = 10;
A_temp = sys.Box_x;
b_temp = sys.box_x;
for i = 1:nSteps
    
    % shrink the disturbance space
    for j = 1:length(sys.box_x)
         
        cvx_begin quiet
            variable w(sys.n,1) 
            variable J
            
            maximize J 
            subject to 
                J == sys.Box_x(j,:)*A_nom^i*w;
                sys.H_w*w <= sys.h_w;                
                
        cvx_end
        
        max_val(j,i) = J;
        clear J w
    end    
 A_temp = [A_temp ; sys.Box_x*A_nom^i];
 b_temp = [b_temp ; sys.box_x-sum(max_val(:,1:i),2)];
end

figure; plotregion(-A_temp,-b_temp);
X0.A = A_temp;
X0.b = b_temp;

% Need to compute vertices of X0


end