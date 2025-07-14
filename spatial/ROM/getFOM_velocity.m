function [V] = getFOM_velocity(R,t,options)
% get FOM velocity based on ROM coefficients

B   = options.rom.B;
Vbc = options.rom.Vbc;

if options.rom.lifting_fun
    V = B*R + options.rom.F_inhom*get_a_bc(t,options);
else
    V = B*R + Vbc;
end
