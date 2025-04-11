function options = diffusiveGersh(options)
  visc = options.case.visc;
  Omu_inv = options.grid.Omu_inv;
  Omv_inv = options.grid.Omv_inv;
  Nu = options.grid.Nu;
  Nv = options.grid.Nv;
  Diffu = options.discretization.Diffu;
  Diffv = options.discretization.Diffv;
  switch visc
    case 'laminar'
    test = spdiags(Omu_inv,0,Nu,Nu)*Diffu;
    sum_diff_u = abs(test)*ones(Nu,1) - diag(abs(test)) - diag(test);
    test = spdiags(Omv_inv,0,Nv,Nv)*Diffv;
    sum_diff_v = abs(test)*ones(Nv,1) - diag(abs(test)) - diag(test);
    eb = max([max(sum_diff_u) max(sum_diff_v)]);
    otherwise
      error('other diffusive methods than laminar not implemented.')
  end
  options.stability.gersh.eb_D = eb;
end
