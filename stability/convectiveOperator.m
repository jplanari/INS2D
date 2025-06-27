function [Ccu, Ccv] = convectiveOperator(C, options, order4)
    indu = options.grid.indu;
    indv = options.grid.indv;
    
    N1 = options.grid.N1;
    N2 = options.grid.N2;
    N3 = options.grid.N3;
    N4 = options.grid.N4;

    uh = C(indu);
    vh = C(indv);

    if (order4 == 0)
        Cux = options.discretization.Cux;
        Cuy = options.discretization.Cuy;
        Cvx = options.discretization.Cvx;
        Cvy = options.discretization.Cvy;
    
        Au_ux = options.discretization.Au_ux;
        Au_uy = options.discretization.Au_uy;
        Av_vx = options.discretization.Av_vx;
        Av_vy = options.discretization.Av_vy;
    
        Iu_ux = options.discretization.Iu_ux;
        Iv_uy = options.discretization.Iv_uy;
        Iu_vx = options.discretization.Iu_vx;
        Iv_vy = options.discretization.Iv_vy;
    
        yIu_ux = options.discretization.yIu_ux;
        yIv_uy = options.discretization.yIv_uy;
        yIu_vx = options.discretization.yIu_vx;
        yIv_vy = options.discretization.yIv_vy;
    else
         error('order4 implementation in evaluateEV.m not done');
    end

    Ccu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
         Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
    Ccv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
         Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy; 
end