function [options] = evaluateEV_ROM(R, options)
    B = options.rom.B;
    C = B*R;

    [Ccu, Ccv] = convectiveOperator(C, options, 0);
    
    Diff = options.rom.Diff;

    Conv = B'*blkdiag(Ccu,Ccv)*B;
    
    ev_C = eig(Conv);
    ev_D = eig(Diff);
    ev_CD = eig(Diff-Conv);

    options.rom.stability.ev_C = ev_C;
    options.rom.stability.ev_D = ev_D;
    options.rom.stability.ev_CD = ev_CD;
    options.rom.stability.eb_C = max(abs(ev_C));
    options.rom.stability.eb_D = max(abs(ev_D));
end

