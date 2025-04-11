%% ANALYZE PERFORMANCE OF GERSHGORIN vs AlgEigCD

clc;
close all;
clear all;

addpath('stability\');
addpath('stability/TG_options/');

cases = [10, 20, 40, 80];

times_AECD_conv = zeros(size(cases));
times_AECD_diff = zeros(size(cases));
times_gers_conv = zeros(size(cases));
times_gers_diff = zeros(size(cases));

nRuns = 10000;

for i=1:length(cases)
    ms = cases(i);
    load(strcat('TG_options_',num2str(ms)));
    t1=0;
    t2=0;
    t3=0;
    t4=0;
    options = setupAECD(options);
    dummyVel = rand(options.grid.NV,1);
    for j=1:nRuns
        tS=tic; options=diffusiveGersh(options); t1 = t1+toc(tS);
        tS=tic; options=diffusiveAECD(options); t2 = t2+toc(tS);
        tS=tic; options=convectiveGersh(dummyVel,options); t3 = t3+toc(tS);
        tS=tic; options=convectiveAECD(dummyVel,options); t4 = t4+toc(tS);
    end

    times_gers_diff(i) = t1/nRuns;
    times_AECD_diff(i) = t2/nRuns;
    times_gers_conv(i) = t3/nRuns;
    times_AECD_conv(i) = t4/nRuns;
end

figure(1)
plot(cases,times_AECD_conv*1e3,'-xk','MarkerSize',15);
hold on
grid on
plot(cases,times_gers_conv*1e3,'-ok','MarkerSize',15);
plot(cases,times_AECD_diff*1e3,'--xk','MarkerSize',15);
plot(cases,times_gers_diff*1e3,'--ok','MarkerSize',15);
xticks(cases);
xlabel('$N$','Interpreter','latex')
ylabel('Wall-clock time [ms]','Interpreter','latex')
legend({'$\rho(C)$, \texttt{AlgEigCD}','$\rho(C)$, Gershgorin', ...
    '$\rho(D)$, \texttt{AlgEigCD}', '$\rho(D)$, Gershgorin'},...
    'Interpreter','latex','Location','best')

figure(2)
plot(cases,times_gers_conv./times_AECD_conv,'-xk','MarkerSize',15);
hold on
grid on
plot(cases,times_gers_diff./times_AECD_diff,'--xk','MarkerSize',15);
xticks(cases);
xlabel('$N$','Interpreter','latex')
ylabel('Speed-up','Interpreter','latex')
legend({'$\rho(C)$','$\rho(D)$'},'Interpreter','latex','Location','best')