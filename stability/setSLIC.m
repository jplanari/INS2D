function [U] = setSLIC(options)
    yu = options.grid.yu;
    xv = options.grid.xv;
    
    d   = pi/15;
    e   = 0.05;
    u   = tanh( (yu-pi/2)/d) .* (yu<=pi) + tanh( (3*pi/2 - yu)/d) .* (yu>pi);
    v   = e*sin(xv);
    U = [u(:);v(:)];
end

