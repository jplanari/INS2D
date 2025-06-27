function [Fsxx,Fsxy,Fsyx,Fsyy] = computeFluxes(V,options)

    indu = options.grid.indu;
    indv = options.grid.indv;

    Iu_ux = options.discretization.Iu_ux;
    Iv_uy = options.discretization.Iv_uy;
    Iu_vx = options.discretization.Iu_vx;
    Iv_vy = options.discretization.Iv_vy;

    yIu_ux = options.discretization.yIu_ux;
    yIv_uy = options.discretization.yIv_uy;
    yIu_vx = options.discretization.yIu_vx;
    yIv_vy = options.discretization.yIv_vy;

    % N1 = options.grid.N1;
    % N2 = options.grid.N2;
    % N3 = options.grid.N3;
    % N4 = options.grid.N4;

    uh = V(indu);
    vh = V(indv);

    % Fsxx = spdiags(Iu_ux*uh+yIu_ux,0,N1,N1);
    % Fsxy = spdiags(Iv_uy*vh+yIv_uy,0,N2,N2);
    % Fsyx = spdiags(Iu_vx*uh+yIu_vx,0,N3,N3);
    % Fsyy = spdiags(Iv_vy*vh+yIv_vy,0,N4,N4);

    Fsxx = Iu_ux*uh+yIu_ux;
    Fsxy = Iv_uy*vh+yIv_uy;
    Fsyx = Iu_vx*uh+yIu_vx;
    Fsyy = Iv_vy*vh+yIv_vy;
    
end

