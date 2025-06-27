function [options] = evaluateCrEV(options,method)
    M = options.rom.M;
    evCR = zeros(M,1);
    Cq = options.rom.Conv_quad;
    for i = 1:M
        cols = i:M:M^2-M+i;
        Ci = Cq(:,cols);
        if strcmp(method,'eig')
            evCR(i) = max(abs(eig(Ci)));
        else
            evCR(i) = gershgorin(Ci);
        end
    end
    options.rom.ev_Cr = evCR;
end

