function options = setupAECD_stg(options)
    order4 = options.discretization.order4;
    isrom = options.rom.rom;
    Nu = options.grid.Nu;
    Nv = options.grid.Nv;

    if(order4==0)
        Cux = options.discretization.Cux;
        Cuy = options.discretization.Cuy;
        Cvx = options.discretization.Cvx;
        Cvy = options.discretization.Cvy;
    
        Su_ux = options.discretization.Su_ux;
        Su_uy = options.discretization.Su_uy;
        Sv_vx = options.discretization.Sv_vx;
        Sv_vy = options.discretization.Sv_vy;

        Au_ux = options.discretization.Au_ux;
        Au_uy = options.discretization.Au_uy;
        Av_vx = options.discretization.Av_vx;
        Av_vy = options.discretization.Av_vy;


    else
        error('fuck!! order4 implementation for AlgEigCD is not done!')
    end
    
    Cu = [Cux Cuy];
    Su = [Su_ux; Su_uy];
    Au = [Au_ux; Au_uy];
    Cv = [Cvx Cvy];
    Sv = [Sv_vx; Sv_vy];
    Av = [Av_vx; Av_vy];
    
    Oinv_u = spdiags(options.grid.Omu_inv,0,Nu,Nu);
    Oinv_v = spdiags(options.grid.Omv_inv,0,Nv,Nv);

    % C = blkdiag(Cu,Cv);
    % A = blkdiag(Au,Av);
    % O = blkdiag(Oinv_u,Oinv_v);

    AECDu = abs(abs(Cu)'*Oinv_u*abs(Au)');
    AECDv = abs(abs(Cv)'*Oinv_v*abs(Av)');

    AECDu_d = abs(Cu'*Oinv_u*Su');
    AECDv_d = abs(Cv'*Oinv_v*Sv');

    options.stability.AECDu = AECDu;
    options.stability.AECDv = AECDv;
    options.stability.AECDu_d = AECDu_d;
    options.stability.AECDv_d = AECDv_d;
    % options.stability.AECD_s = abs(C'*O*A');
    % options.stability.AECD_s = abs(C')*O*A';


end

