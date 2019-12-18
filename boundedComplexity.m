function cont = boundedComplexity(cont)
%% add constraints on parameters with bounded complexity overall

%% Decide on directions to add constraints
% number of parameters
p = size(cont.H_theta,2);

% Each plane in parameter space can be specified by a direction vector
% chosen. The directions involving two zeros and 1 are chosen already. Now,
% we vary each parameter in the space between -a and b, with n_values
% number of grid points in between. 

n_values = 3; 
values = [linspace(-1.1,1,n_values)'];
l_values = length(values); % length of the values vector

new_dir = kron(values,ones(l_values,1));
for i =1:p-1
    new_dir = [kron(new_dir,ones(l_values,1)),kron(ones(l_values^(i+1),1),values)];
end

n_new_dir = size(new_dir,1);
%% Decide on lower and upper bounds for these directions

% Check if already calculated bounds are valid
try
    % load and check length
    matrices = load('new_bounds.mat');
    calc_flag = length(matrices.lb_new)~=n_new_dir;
    
    % check equality
    if any(any(new_dir-matrices.new_dir))
        % not same matrices
        calc_flag = 1;
    end
    
catch
    calc_flag = 1;
end


if calc_flag
    % calculate upper and lower bounds again
    lb_new = NaN*ones(n_new_dir,1);
    ub_new = NaN*ones(n_new_dir,1);
    for i = 1:n_new_dir 
    % calculate upper bound   
        cvx_begin quiet
        cvx_solver gurobi
            variable x_max(p,1)
            variable J_max
            maximize J_max

            subject to
            J_max == new_dir(i,:)*x_max;
            cont.H_theta*x_max <= cont.h_theta_k;
        cvx_end

        if (isnan(J_max))
           warning(cvx_status) 
           flag = 1;
           return;
        end
        ub_new(i) = J_max;

    % calculate lower bound   
        cvx_begin quiet
        cvx_solver gurobi
            variable x_min(p,1)
            variable J_min
            minimize J_min

            subject to
            J_min == new_dir(i,:)*x_min;
            cont.H_theta*x_min <= cont.h_theta_k;
        cvx_end

        if (isnan(J_min))
           warning(cvx_status) 
           flag = 1;
           return;
        end
        lb_new(i) = J_min;

    end
else % use loaded matrices
    warning('Bounded complexity matrices are being loaded from new_bounds.mat.')
    ub_new = matrices.ub_new;
    lb_new = matrices.lb_new;
end

save('new_bounds.mat','new_dir','lb_new','ub_new');

cont.H_theta = [cont.H_theta; new_dir; -new_dir];
cont.h_theta_k = [cont.h_theta_k; ub_new; -lb_new];
end
