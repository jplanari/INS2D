function options = evaluateEV_FOM(C, options)
    order4 = options.discretization.order4;
    Ouinv = sparse(diag(options.grid.Omu_inv));
    Ovinv = sparse(diag(options.grid.Omv_inv));
    [Ccu, Ccv] = convectiveOperator(C, options, order4);
    [Dcu, Dcv] = diffusiveOperator(options);
    
    ev_C = [eig(full(Ouinv*Ccu)); eig(full(Ovinv*Ccv))];
    ev_D = [eig(full(Ouinv*Dcu)); eig(full(Ovinv*Dcv))];
    ev_CD = [eig(full(Ouinv*(-Ccu+Dcu))); eig(full(Ovinv*(-Ccv+Dcv)))];

    options.stability.ev_C = ev_C;
    options.stability.ev_D = ev_D;
    options.stability.ev_CD = ev_CD;
    options.stability.eb_C = max(imag(ev_CD));
    options.stability.eb_D = max(abs(real(ev_D)));
end



