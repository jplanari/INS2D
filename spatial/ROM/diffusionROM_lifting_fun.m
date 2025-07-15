function [d2, Jac] = diffusionROM_lifting_fun(R,t,options,getJacobian)

visc = options.case.visc;

Diff  = options.rom.Diff;

switch visc

    case 'laminar'
        d2     = Diff*R;

        R_bc    = get_a_bc(t,options);
        yDiff = options.rom.yDiff;

        d2 = d2 + yDiff*R_bc;

        Jac    = Diff;
    otherwise
        error('turbulent diffusion not implemented for ROM');

end

end

