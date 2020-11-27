function D = D_mult(sys,x,u)
    x_mult = sys.Ap_mult*x;
    x_mult = reshape(x_mult,sys.n,[]);
    u_mult = sys.Bp_mult*u;
    u_mult = reshape(u_mult,sys.n,[]);    
    D = x_mult + u_mult;
end
