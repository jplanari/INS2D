function F_inhom = get_F_inhom(phi_bc,options)

F_M = options.discretization.F_M;
Y_M = F_M*phi_bc;

L = options.discretization.A;
Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
G = [Gx;Gy];
Om_inv = options.grid.Om_inv;
F_inhom = -Om_inv.*(G*(L\Y_M)); %pfusch