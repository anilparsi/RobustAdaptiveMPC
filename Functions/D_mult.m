function D = D_mult(sys,x,u)
    x_mult = reshape(sys.Ap,[sys.n*sys.p,sys.n])*x;
    x_mult = reshape(x_mult,[sys.n,sys.p]);
    u_mult = reshape(sys.Bp,[sys.n*sys.p,sys.m])*u;
    u_mult = reshape(u_mult,[sys.n,sys.p]);    
    D = x_mult + u_mult;
end