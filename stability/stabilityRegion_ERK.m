function [options] = stabilityRegion_ERK(options)

r = linspace(0.1,5,10000);
phi = linspace(0,pi/2,1000);
rStab = zeros(size(phi));

[~,b,~,~] = getRKmethod(options.time.RK);
s = length(b);

for i=1:length(phi)
    for j=1:length(r)
        zz = -r(j)*cos(phi(i))+1i*r(j)*sin(phi(i));
        R = stabFun(s,zz);
        if(abs(abs(R)-1.0) < 1e-3)
            rStab(i)=r(j);
            break;
        end
    end
end

options.stability.region.phi=phi;
options.stability.region.R=rStab;

end


