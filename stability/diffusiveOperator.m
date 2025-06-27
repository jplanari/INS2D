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