function [yMt] = determine_yMt(tn,options)
    t_snp = options.rom.t_snp;
    yMt_full = options.rom.yMt;
    idx_1 = find(t_snp > tn, 1, 'first');
    idx_2 = find(t_snp < tn, 1, 'last');
    if (isempty(idx_2))
        idx_2 = 1;
    end
    if (isempty(idx_1))
        idx_1 = length(t_snp);
    end

    yMt = (t_snp(idx_2)-tn)/(t_snp(idx_2)-t_snp(idx_1))*yMt_full(:,idx_1)...
        + (-t_snp(idx_1)+tn)/(t_snp(idx_2)-t_snp(idx_1))*yMt_full(:,idx_2);

end

