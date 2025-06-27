function [options] = diffusiveAECD(options)
    visc = options.case.visc;
    AECD = options.stability.AECD;
    switch visc
        case 'laminar'
            nu = 1/options.fluid.Re;
        otherwise
            error('other diffusive methods than laminar not implemented yet')
    end
    Ds = options.discretization.Ds;
    vec = nu*Ds;

    eb = max(AECD*vec);

    options.stability.bounds.eb_D = eb;
end

