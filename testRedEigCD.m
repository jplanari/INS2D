%% METHOD QUALITY ASSESSMENT FOR RedEigCD

clc;
close all;
clear;

addpath('stability\');
addpath('stability/SL_ROM_options/');


M = [2,4,8,16,32];

eigm_Dr = zeros(size(M));
gersh_Dr = zeros(size(M));
eigm_Cr = zeros(size(M));
gersh_Cr = zeros(size(M));
eigma_Cr = zeros(size(M));
gersha_Cr = zeros(size(M));
t_e_Dr = zeros(size(M));
t_g_Dr = zeros(size(M));
t_e_Cr = zeros(size(M));
t_g_Cr = zeros(size(M));
t_ea_Cr = zeros(size(M));
t_ga_Cr = zeros(size(M));
nRuns = 10000;

for i=1:length(M)
    m = M(i);
    load(strcat('SL_ROM_options_',num2str(m)));
    B = options.rom.B;
    V = setSLIC(options);
    a = B'*V;
    Dr = options.rom.Diff;
    Cq = options.rom.Conv_quad;
    for k = 1:nRuns
    %eigm(Dr)
    tS = tic; eigm_Dr(i) = max(abs(eig(Dr))); t_e_Dr(i) = t_e_Dr(i)+toc(tS);
    %gersh(Dr)
    tS = tic; gersh_Dr(i) = gershgorin(Dr); t_g_Dr(i) = t_g_Dr(i)+toc(tS);
    %Build a_i*C_ri
    tS = tic;
    C_tens = splitCquad(Cq,m);
    Cl = zeros(m,m);
    for j=1:m
        Cl = Cl + a(j)*C_tens(:,:,j);
    end
    tsplit = toc(tS);
    %eigm(Cr)
    tS = tic; eigm_Cr(i) = max(abs(eig(Cl))); t_e_Cr(i) = t_e_Cr(i)+toc(tS);
    %gersh(Cr)
    tS = tic; gersh_Cr(i) = gershgorin(Cl); t_g_Cr(i) = t_g_Cr(i)+toc(tS);
    %eigma(Cr)
    options = evaluateCrEV(options,'eig');
    tS = tic; eigma_Cr(i) = abs(a)'*options.rom.ev_Cr; t_ea_Cr(i) = t_ea_Cr(i)+toc(tS);
    %gersha(Cr)
     options = evaluateCrEV(options,'gersh');
    tS = tic; gersha_Cr(i) = abs(a)'*options.rom.ev_Cr; t_ga_Cr(i) = t_ga_Cr(i)+toc(tS);              
    end
end

e_g_Dr = abs(eigm_Dr-gersh_Dr)./eigm_Dr;
e_g_Cr = abs(eigm_Cr-gersh_Cr)./eigm_Cr;
e_ea_Cr = abs(eigm_Cr-eigma_Cr)./eigm_Cr;
e_ga_Cr = abs(eigm_Cr-gersha_Cr)./eigm_Cr;

%Diffusive eigenbounds
figure;
semilogx(M,eigm_Dr,'-xk','MarkerSize',15)
hold on
grid on
semilogx(M,gersh_Dr,'-ok','MarkerSize',15)
xticks(M);
xlim([M(1),M(end)])
xlabel('$M$','Interpreter','latex');
ylabel('$\rho(D_r)$','Interpreter','latex');
legend({'\texttt{eigm}($D_r$)','\texttt{gersh}($D_r$)'},...
    'interpreter','latex');

%Convective eigenbounds
figure;
semilogx(M,eigm_Cr,'-xk','MarkerSize',15)
hold on
grid on
semilogx(M,gersh_Cr,'-ok','MarkerSize',15)
semilogx(M,eigma_Cr,'-sk','MarkerSize',15)
semilogx(M,gersha_Cr,'-dk','MarkerSize',15)
xticks(M);
xlim([M(1),M(end)])
xlabel('$M$','Interpreter','latex');
ylabel('$\rho(C_r)$','Interpreter','latex');
legend({'\texttt{eigm}($C_r$)','\texttt{gersh}($C_r$)',...
    '\texttt{eigma}($C_r$)','\texttt{gersha}($C_r$)'},...
    'interpreter','latex');

%Wall-clock time, diffusive
figure;
loglog(M,t_e_Dr/nRuns*1000,'-xk','MarkerSize',15)
hold on
grid on
loglog(M,t_g_Dr/nRuns*1000,'-ok','MarkerSize',15)
xticks(M);
xlim([M(1),M(end)])
xlabel('$M$','Interpreter','latex');
ylabel('$t$ [ms]','Interpreter','latex');
legend({'\texttt{eigm}($D_r$)','\texttt{gersh}($D_r$)'},...
    'interpreter','latex');


%Wall-clock time, convective
figure;
loglog(M,t_e_Cr/nRuns*1000,'-xk','MarkerSize',15)
hold on
grid on
loglog(M,t_g_Cr/nRuns*1000,'-ok','MarkerSize',15)
loglog(M,t_ea_Cr/nRuns*1000,'-sk','MarkerSize',15)
loglog(M,t_ga_Cr/nRuns*1000,'-dk','MarkerSize',15)
xticks(M);
xlim([M(1),M(end)])
xlabel('$M$','Interpreter','latex');
ylabel('$t$ [ms]','Interpreter','latex');
legend({'\texttt{eigm}($C_r$)','\texttt{gersh}($C_r$)',...
    '\texttt{eigma}($C_r$)','\texttt{gersha}($C_r$)'},...
    'interpreter','latex');

%Method errors vs exact spectral radius
figure;
loglog(M,e_g_Dr,'-xk','MarkerSize',15)
hold on
grid on
loglog(M,e_g_Cr,'-ok','MarkerSize',15)
loglog(M,e_ea_Cr,'-sk','MarkerSize',15)
loglog(M,e_ga_Cr,'-dk','MarkerSize',15)
xticks(M);
xlim([M(1),M(end)])
ylim([1e-5 10])
xlabel('$M$','Interpreter','latex');
ylabel('$\varepsilon$','Interpreter','latex');
legend({'\texttt{gersh}($D_r$)','\texttt{gersh}($C_r$)',...
    '\texttt{eigma}($C_r$)','\texttt{gersha}($C_r$)'},...
    'interpreter','latex');
