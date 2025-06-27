function [options] = setupAECD(options)
%Incidence matrices
G = options.discretization.G;
Tcs = G./abs(G);
Tcs(isnan(Tcs))=0;

%Collocated grid volume inverse, Identity for ROM
Np = options.grid.Np;

if(options.rom.rom==1)
    B = options.rom.B;
    Oinv = speye(Np);
    Tcs = B'*Tcs;

else
    Omp_inv = options.grid.Omp_inv;
    Oinv = spdiags(Omp_inv,0,Np,Np);
end

Tsc = Tcs';

%Vector of areas
hxi = options.grid.hxi;
hyi = options.grid.hyi;
Nux_in = options.grid.Nux_in;
Nvy_in = options.grid.Nvy_in;

As = [repmat(hyi(:), Nux_in, 1); repmat(hxi(:), Nvy_in, 1)];
Ds = [repmat(hxi(:), Nux_in, 1); repmat(hyi(:), Nvy_in, 1)];

%Relevant matrix
AECD = abs(Tcs*Oinv*Tsc);

%Assign to options struct
options.discretization.Tsc = Tsc;
options.discretization.Tcs = Tcs;
options.discretization.As = As;
options.discretization.Ds = As./Ds(1:length(As));
options.stability.AECD = AECD;


end

