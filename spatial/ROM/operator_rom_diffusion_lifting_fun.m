function [yDiff,Diff] = operator_rom_diffusion_lifting_fun(P,options)
% yDiff: linear in a_bc
% Diff: linear in a_hom

B   = options.rom.B;

Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;

Diff_   = P*blkdiag(Diffu, Diffv)*B;

phi_bc = options.rom.phi_bc;
M_bc = size(phi_bc,2);
M = options.rom.M;
yDiff_BC_ = zeros(M,M_bc);

for i = 1:M_bc
    yBC = phi_bc(:,i);
    options = set_bc_vectors_from_yBC(options,yBC);

    yDiffu1 = options.discretization.yDiffu;
    yDiffv1 = options.discretization.yDiffv;
    yDiff_BC_(:,i) = P*[yDiffu1; yDiffv1];
end


F_inhom = options.rom.F_inhom;

Diff_inhom_ = P*blkdiag(Diffu, Diffv)*F_inhom;



Diff   = Diff_;
yDiff = yDiff_BC_ + Diff_inhom_;



