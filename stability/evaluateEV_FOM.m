function options = evaluateEV_FOM(C, options)
    order4 = options.discretization.order4;
    Oinv = sparse(diag(options.grid.Omp_inv));
    [Ccu, Ccv] = convectiveOperator(C, options, order4);
    [Dcu, Dcv] = diffusiveOperator(options);
    
    ev_C = [eig(full(Oinv*Ccu)); eig(full(Oinv*Ccv))];
    ev_D = [eig(full(Oinv*Dcu)); eig(full(Oinv*Dcv))];
    ev_CD = [eig(full(Oinv*(-Ccu+Dcu))); eig(full(Oinv*(-Ccv+Dcv)))];

    options.stability.ev_C = ev_C;
    options.stability.ev_D = ev_D;
    options.stability.ev_CD = ev_CD;
    options.stability.eb_C = max(abs(ev_C));
    options.stability.eb_D = max(abs(ev_D));
end

function [Dcu, Dcv] = diffusiveOperator(options)
    visc = options.case.visc;
    
    switch visc
        case 'laminar'
            Dcu = options.discretization.Diffu;
            Dcv = options.discretization.Diffv;
        otherwise
            error('other diffusive methods than laminar not implemented')
    end
end

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