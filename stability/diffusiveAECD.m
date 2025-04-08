function options = diffusiveAECD(options)
   visc = options.case.visc;
   AECD = options.discretization.AECD;
   switch visc
       case 'laminar'
           diagLs = ones(size(AECD,1),1)* (1/options.fluid.Re);
           eb = max(abs(AECD)*diagLs);
       otherwise
           error('other diffusive methods than laminar not implemented yet.')
   end
   options.stability.AECD.eb_D = eb;
end