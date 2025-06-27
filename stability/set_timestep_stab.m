function [dt,options] = set_timestep_stab(V,options)

    %% Estimate eigenbounds with AlgEigCD. Requires options=setupAECD(options) before the time loop.

    [Fsxx,Fsxy,Fsyx,Fsyy] = computeFluxes(V,options);
    Fsx = abs([Fsxx; Fsxy]);
    Fsy = abs([Fsyx; Fsyy]);
    options = convectiveAECD_stg(Fsx,Fsy,options);

    % options = diffusiveGersh(options);
    % options = convectiveGersh(V,options);

    [~,b,~,~] = getRKmethod(options.time.RK);
    s = length(b);
    
    ebR = options.stability.bounds.eb_D;
    ebI = options.stability.bounds.eb_C;

    eb_n = sqrt(ebR^2+ebI^2);
    phi = atan(ebI/ebR);
    
    r = linspace(0.1,5,10000);
    for i=1:length(r)
        zz = -r(i)*cos(phi)+1i*r(i)*sin(phi);
        R = stabFun(s,zz);
        if(abs(abs(R)-1.0) < 1e-3)
            rStab=r(i);
            break;
        end
    end

    dt = rStab/eb_n;
end