function options = convectiveGersh(C, options)
  order4 = options.discretization.order4;
  indu = options.grid.indu;
  indv = options.grid.indv;

  N1 = options.grid.N1;
  N2 = options.grid.N2;
  N3 = options.grid.N3;
  N4 = options.grid.N4;
  Nu = options.grid.Nu;
  Nv = options.grid.Nv;

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
    
    Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
         Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
    Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
         Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy; 

    test = spdiags(options.grid.Omu_inv,0,Nu,Nu)*Cu;
    %sum_conv_u = abs(test)*ones(Nu,1) - diag(abs(test)) - diag(test);
    sum_conv_u = abs(test)*ones(Nu,1);
    test = spdiags(options.grid.Omv_inv,0,Nv,Nv)*Cv;
    %sum_conv_v = abs(test)*ones(Nv,1) - diag(abs(test)) - diag(test);
    sum_conv_v = abs(test)*ones(Nv,1);
    options.stability.gersh.eb_C = max([max(sum_conv_u) max(sum_conv_v)]);

end
