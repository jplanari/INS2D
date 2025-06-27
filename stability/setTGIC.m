function [U] = setTGIC(options)
    xu = options.grid.xu;
    yu = options.grid.yu;
    xv = options.grid.xv;
    yv = options.grid.yv;

    u  = - sin(pi*xu).*cos(pi*yu);
    v  = cos(pi*xv).*sin(pi*yv);
    U = [u(:);v(:)];
end

