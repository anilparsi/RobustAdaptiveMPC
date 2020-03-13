function cont = cont_params(sys)
% prediction horizon
cont.N = 8;

% cost matrices
cont.Q = eye(sys.n);
cont.R = eye(sys.m);
cont.P = [2.8  0.56;
          0.56 2.5];
cont.Q_L = chol(cont.Q);
cont.R_L = chol(cont.R);
      
      
% block size for updating h_theta_k
cont.blk = 10;

% Define state tube shape : H_x*x <= vec_1_x, vertices: x_v(:,i)
cont.H_x = [eye(sys.n); -eye(sys.n)];
cont.vec_1_x = ones(length(cont.H_x),1);
cont.nHx = size(cont.H_x,1);  % denoted u in Lorenzen(2019)
cont.x_v = [1 1; 1 -1; -1 1; -1 -1 ]';
cont.nx_v = length(cont.x_v); % number of vertices

% prestabilizing gain
cont.w_bar = compute_wbar(sys,cont);
[cont.K,cont.alpha_bar,cont.alpha_min] = prestab_controller(sys,cont);

% Define terminal constraint: z_N|k = 0; h_T*alpha_N|k <= 1;
cont.f_bar = max((sys.F+sys.G*cont.K)*cont.x_v,[],2);

% Bounded complexity
cont.H_theta = sys.H_theta;
cont.h_theta_k = sys.h_theta;
cont.h_theta_0 = sys.h_theta;

end