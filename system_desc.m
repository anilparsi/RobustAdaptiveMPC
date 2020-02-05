function sys = system_desc()

% define matrices
sys.A0 = [ 0.5 0.2;
      -0.1 0.6];
sys.Ap(:,:,1) = [0.042 0
             0.072 0.03];
sys.Ap(:,:,2) = [0.015 0.019
             0.009 0.035];
sys.Ap(:,:,3) = zeros(2,2);

sys.B0 = [0; 0.5];
sys.Bp(:,:,1) = zeros(2,1);
sys.Bp(:,:,2) = zeros(2,1);
sys.Bp(:,:,3) = [0.04; 0.054];

% define dimensions
sys.n = size(sys.Bp,1);
sys.m = size(sys.Bp,2);
sys.p = size(sys.Bp,3);

% define bounds on theta: H_theta*theta <= h_theta
sys.H_theta = [eye(sys.p);-eye(sys.p)];
sys.h_theta = [ones(sys.p,1);ones(sys.p,1)];

% Generate H_theta with a predetermined complexity
sys = boundedComplexity(sys);
sys.nHtheta = length(sys.h_theta);

% define disturbance bounds: H_w*w <= h_w
sys.w_bound = 0.1;
sys.H_w = [eye(sys.n);-eye(sys.n)]/sys.w_bound;
sys.h_w = [ones(sys.n,1);ones(sys.n,1)];
sys.nHw = length(sys.h_w);


% define state and input constraints: F*x + G*u <= vec_1_cons
sys.F = [0 -1/0.3; 
         0  0
         0  0];
sys.G = [0; 1; -1];
sys.nc = size(sys.F,1);
sys.vec_1_cons = ones(sys.nc,1);

% define box constraint on x to find mu: Box_x*x <= box_x
sys.Box_x = [eye(sys.n),-eye(sys.n)];
sys.box_x = 3*[ones(sys.n,1);-ones(sys.n,1)];

sys.x0 = [2;3];
end

