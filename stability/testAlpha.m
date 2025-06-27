clear
clc
close all;

addpath('stability\');
addpath('stability/LDC_options/');

alpha = linspace(-0.5,1.5);
% cases = [10,20,40];
cases = [10,20  ,40,80];


linestyles={'-k','--k','-.k',':k'};

nRuns=4;

for i=1:length(cases)
    ms = cases(i);
    load(strcat('LDC_options_',num2str(ms)));
    options = setupAECD_stg(options);
    ebc = zeros(size(alpha));
    ebd = zeros(size(alpha));
    dummyVel = rand(options.grid.NV,1);
    b = options.discretization.M*dummyVel;
    L = options.discretization.M*options.discretization.G;
    p = L\b;
    dummyVel = dummyVel-options.discretization.G*p;
    [Fsxx,Fsxy,Fsyx,Fsyy] = computeFluxes(dummyVel,options);
    Fsx = abs([Fsxx;Fsxy]);
    Fsy = abs([Fsyx;Fsyy]);
    for k=1:nRuns

        for j=  1:length(alpha)
            options.stability.alpha=alpha(j);
            options = convectiveAECD_stg(Fsx,Fsy,options); 
            options = diffusiveAECD_stg(options);
            ebc(j) = ebc(j) + options.stability.bounds.eb_C;
            ebd(j) = ebd(j) + options.stability.bounds.eb_D;
        end
    end
    ebc = ebc/nRuns;
    ebd = ebd/nRuns;
    figure(1);
    semilogy(alpha,ebc,linestyles{i},'MarkerSize',15);
    hold on;
    figure(2);
    semilogy(alpha,ebd,linestyles{i},'MarkerSize',15);
    hold on;
end
figure(1)
grid on
xlabel('$\alpha$','interpreter','latex');
ylabel('$\rho(C)(\alpha)$','interpreter','latex');
legend({'$N=10$','$N=20$','$N=40$','$N=80$'},'interpreter','latex')

figure(2)
grid on
xlabel('$\alpha$','interpreter','latex');
ylabel('$\rho(D)(\alpha)$','interpreter','latex');
legend({'$N=10$','$N=20$','$N=40$','$N=80$'},'interpreter','latex')