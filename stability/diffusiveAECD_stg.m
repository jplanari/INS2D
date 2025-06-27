function [options] = diffusiveAECD_stg(options)
    AECDu = options.stability.AECDu_d;
    AECDv = options.stability.AECDv_d;
    alpha = options.stability.alpha;
    visc = options.case.visc;
    switch visc
        case 'laminar'
            nu = 1/options.fluid.Re;
        otherwise
            error('other diffusive methods than laminar not implemented yet')
    end

    % F = [Fsx; Fsy];

    % nu_u = nu*ones(size(AECDu,1),1);
    % nu_v = nu*ones(size(AECDv,1),1);

    ebu = 0.25*max(nu.^alpha.*AECDu*nu.^(1-alpha));
    ebv = 0.25*max(nu.^alpha.*AECDv*nu.^(1-alpha));
    options.stability.bounds.eb_D = max([ebu,ebv]);
end