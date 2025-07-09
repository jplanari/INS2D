ebI = abs(R')*options.rom.ev_Cr+options.rom.eb_Cl;
% ebI = abs(R')*options.rom.ev_Cr;
ebR = options.rom.eb_D;

[~,b,~,~] = getRKmethod(options.time.RK);
s = length(b);

eb_n = sqrt(ebR^2+ebI^2);
phi = atan(ebI/ebR);

r = linspace(0.1,5,100000);
for i=1:length(r)
    zz = -r(i)*cos(phi)+1i*r(i)*sin(phi);
    R = stabFun(s,zz);
    if(abs(abs(R)-1.0) < 1e-4)
        rStab=r(i);
        break;
    end
end

dt = rStab/eb_n;