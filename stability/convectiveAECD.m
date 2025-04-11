function [options] = convectiveAECD(V,options)
    AECD = options.stability.AECD;
    As = options.discretization.As;
    Fs = abs(V.*As);
    eb = 0.25*max(AECD*Fs);
    options.stability.AlgEigCD.eb_C = eb;
end

