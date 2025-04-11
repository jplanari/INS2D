function [dt,options] = set_timestep_stab(V,options)

    options = diffusiveGersh(options);
    options = convectiveGersh(V,options);
    [~,b,~,~] = getRKmethod(options.time.RK);
    s = length(b);

    eb_n = sqrt(options.stability.gersh.eb_C^2+options.stability.gersh.eb_D^2);
    phi = atan(options.stability.gersh.eb_C/options.stability.gersh.eb_D);
    
    r = linspace(0.1,5,10000);
    for i=1:length(r)
        zz = -r(i)*cos(phi)+1i*r(i)*sin(phi);
        R = stabFun(s,zz);
        if(abs(abs(R)-1.0) < 1e-3)
            rStab=r(i);
            break;
        end
    end

    dt = rStab/eb_n;
end