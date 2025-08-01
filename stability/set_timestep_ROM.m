Clt = options.rom.Conv_linear*kron(eye(M),get_a_bc(t,options));
ebI = abs(R')*options.rom.ev_Cr+gershgorin(0.5*(Clt-Clt'));
ebR = options.rom.eb_D+gershgorin(0.5*(Clt+Clt'));

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
