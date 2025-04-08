function options = diffusiveGersh(options)
  visc = options.case.visc;
  Ouinv = sparse(diag(options.grid.Omu_inv));
  Ovinv = sparse(diag(options.grid.Omv_inv));
  switch visc
    case 'laminar'
      Dcu = Ouinv*options.discretization.Diffu;
      Dcv = Ovinv*options.discretization.Diffv;
      eb = max(gershgorin(Dcu),gershgorin(Dcv));
    otherwise
      error('other diffusive methods than laminar not implemented.')
  end
  options.stability.gersh.eb_D = eb;
end
