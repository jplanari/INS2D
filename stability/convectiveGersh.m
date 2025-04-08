function options = convectiveGersh(C, options)
  order4 = options.discretization.order4;
  indu = options.grid.indu;
  indv = options.grid.indv;

  cu = C(indu);
  cv = C(indv);

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

  uf_ux = Iu_ux*cu+yIu_ux;
  vf_uy  = Iv_uy*cv+yIv_uy;
  uf_vx  = Iu_vx*cu+yIu_vx;
  vf_vy  = Iv_vy*cv+yIv_vy;

  Ouinv = sparse(diag(options.grid.Omu_inv));
  Ovinv = sparse(diag(options.grid.Omv_inv));
  
  Ccu = Ouinv*(Cux*sparse(diag(uf_ux))*Au_ux + Cuy*sparse(diag(vf_uy))*Au_uy);
  Ccv = Ovinv*(Cvx*sparse(diag(uf_vx))*Av_vx + Cvy*sparse(diag(vf_vy))*Av_vy);
  
  options.stability.gersh.eb_C = max(gershgorin(Ccu),gershgorin(Ccv));

end
