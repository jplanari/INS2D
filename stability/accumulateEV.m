function [conv,diff] = accumulateEV(conv,diff,options)
    conv = conv + options.stability.bounds.eb_C;
    diff = diff + options.stability.bounds.eb_D;
end

