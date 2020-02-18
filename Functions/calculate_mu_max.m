function mu = calculate_mu_max(sys)
% 1/mu >= sup_(x,u)  ||D(x,u)||^2
% uses the corner points of x and u to get the supremum (approximation).
% The optimization to be solved is concave. Hence this apx is used.
% This gives the maximum value of mu that can be used for the system

x_test = [];
u_test = [];
D_val  = [];
for i = 1:length(sys.Box_x_v)
    for j = 1:length(sys.Box_u_v)
        x_test = [x_test sys.Box_x_v(:,i)];
        u_test = [u_test sys.Box_u_v(:,j)];
        D_val(i,j) = norm(D_mult(sys,  sys.Box_x_v(:,i), sys.Box_u_v(:,j)),2)^2;
    end
end
D_max = max(max(D_val));

% mu <= 1/D_max
mu = 1/D_max;
end