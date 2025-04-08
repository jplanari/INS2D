function options = setup_AlgEigCD_eigenbounds(options)
    G = options.discretization.Gx;
    Tcs = G./abs(G);
    Tcs(isnan(Tcs))=0;
    if (options.rom.rom==0)
        Omv_inv = sparse(diag(options.grid.Omv_inv));
        %Omv_inv = sparse(eye(size(Tcs,1)));
    else
        Omv_inv = sparse(eye(size(Tcs,1)));
    end
    options.discretization.AECDu = abs(2*Tcs'*Omv_inv*Tcs); %2* because we need to consider the contribution of both x and y faces
end