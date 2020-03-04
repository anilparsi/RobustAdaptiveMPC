function w_bar =  compute_wbar(sys,cont)
%% Computation of w_bar
nHw = size(sys.H_w,1);
nHx = size(cont.H_x,1);
% max Hx_i w, s.t. H_w w <= h_w
% min h_w' li, s.t. H_w' li == Hx_i, li>=0
lam = sdpvar(nHw,nHx,'full');
Constraints = [lam>=0];
J = 0;
for i = 1:nHx
    Constraints = [Constraints,...
        sys.H_w'*lam(:,i)==cont.H_x(i,:)'];
    J = J+sys.h_w'*lam(:,i);
end
options = sdpsettings('solver','gurobi','verbose',0);
diagnostics = optimize(Constraints,J,options);
if diagnostics.problem
   error(['Calculating w_bar:',diagnostics.info]);
end

for i = 1:nHx
    w_bar(i,1) = value(sys.h_w'*lam(:,i));
end


end