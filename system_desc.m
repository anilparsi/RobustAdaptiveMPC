function sys = system_desc()

% define matrices
sys.A0 = [ 0.85  0.5;
           0.6  0.6 ];
sys.Ap(:,:,1) = [0.1 0
                 0    0.3];
sys.Ap(:,:,2) = [-0. 0.2
                  0.0 0.];

sys.B0 = [1.5;
          0.4];
sys.Bp(:,:,1) = zeros(2,1);
sys.Bp(:,:,2) = zeros(2,1);

% define dimensions
sys.n = size(sys.Bp,1);
sys.m = size(sys.Bp,2);
sys.p = size(sys.Bp,3);

% define bounds on theta: H_theta*theta <= h_theta
sys.H_theta = [eye(sys.p);-eye(sys.p)];
sys.h_theta = [ones(sys.p,1);ones(sys.p,1)];
sys.H_theta_v = [1 1; 1 -1; -1 1;  -1 -1]';


% Generate H_theta with a predetermined complexity
sys = boundedComplexity(sys);
sys.nHtheta = length(sys.h_theta);

% define disturbance bounds: H_w*w <= h_w
sys.w_bound = 0.05;
sys.H_w = [eye(sys.n);-eye(sys.n)]/sys.w_bound;
sys.h_w = [ones(sys.n,1);ones(sys.n,1)];
sys.nHw = length(sys.h_w);


% define state and input constraints: F*x + G*u <= vec_1_cons
sys.F = [-0.1 0;
         0 -0.1; 
         0  0
         0  0];
sys.G = [0; 0; .5; -1/2];
sys.nc = size(sys.F,1);
sys.vec_1_cons = ones(sys.nc,1);

% define box constraint on x to find mu: Box_x_v: vertices of box
sys.Box_x = [eye(sys.n);-eye(sys.n)];
sys.box_x = 5*[ones(sys.n,1);ones(sys.n,1)];
sys.Box_x_v = 5*[1,1; 1,-1; -1,1; -1,-1]';
sys.Box_u_v = [1; -1]';

sys.x0 = [1;1.5];
end

