function [options] = compareJacobians(R, tn, options)

M = options.rom.M;
CR = options.rom.Cr;

if options.rom.rom_bc == 2
    Clt = options.rom.Conv_linear*kron(eye(M),get_a_bc(tn,options));
else
    Clt = zeros(M,M);
end

% Clt = zeros(M,M);

symClt = 0.5*(Clt-Clt');

Jp = symClt;
J = symClt;
for i=1:M
    J(:,i) = CR(:,:,i)*R; %Dropped term
    Jp = Jp + R(i) * CR(:,:,i); %Reduced Jacobian approximation
end

% display(issymmetric(J,"skew"));
% display(issymmetric(Jp,"skew"));

J = J + Jp; %Dropped + reduced

% eig_J = max(abs(norm(eig(J))));
% eig_Jp = max(abs(norm(eig(Jp))));

eig_J = max(abs(imag(eig(J))));
eig_Jp = max(abs(imag(eig(Jp))));


ebI = abs(R')*options.rom.ev_Cr+gershgorin(symClt);
gershI = gershgorin(Jp);

fprintf("%d\t%e\t%e\t%e\t%e\n",M,eig_J,eig_Jp,ebI,gershI);
end