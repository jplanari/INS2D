function [Conv, Jac] = convectionROM_lifting_fun(R,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

if getJacobian    
    M   = options.rom.M;
    E   = speye(M);

    Jac = options.rom.Conv_quad*(kron(E,R)+kron(R,E));
else
    Jac = -666;
end

Conv = options.rom.Conv_quad*kron(R,R);
    
R_bc    = get_a_bc(t,options);

Conv = Conv ...
    + options.rom.Conv_linear*kron(R,R_bc) ...
    + options.rom.yConv*kron(R_bc,R_bc);

if getJacobian
    Jac = Jac ...
        + options.rom.Conv_linear*kron(E,R_bc);
end


