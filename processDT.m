clear
clc
close all

addpath('stability/dt_comparison/full');

M = [16,32,64,128,256];

for i=1:length(M)
    rom_data.(strcat('rom_',num2str(M(i)))) = load(strcat('rom_',strcat(num2str(M(i)),'.mat')));
    innerFields = fieldnames(rom_data.(strcat('rom_',num2str(M(i)))).(strcat('timestep_handle_ROM70_',num2str(M(i)))));

    for j=1:length(innerFields)
        rom_data.(strcat('rom_',num2str(M(i)))).(innerFields{j})=...
            rom_data.(strcat('rom_',num2str(M(i)))).(strcat('timestep_handle_ROM70_',num2str(M(i)))).(innerFields{j});
    end
    rom_data.(strcat('rom_',num2str(M(i)))) = rmfield(rom_data.(strcat('rom_',num2str(M(i)))),strcat('timestep_handle_ROM70_',num2str(M(i))));
end

fom_data=load('fom.mat');
innerFields = fieldnames(fom_data.timestep_handle_FOM70);
for j=1:length(innerFields)
        fom_data.(innerFields{j})=fom_data.timestep_handle_FOM70.(innerFields{j});
end
fom_data = rmfield(fom_data,'timestep_handle_FOM70');

t_fom = zeros(fom_data.n,1);
for i=2:fom_data.n
    t_fom(i) = t_fom(i-1)+fom_data.dts(i);
end

clear innerFields
linestyle = {'--k','-.k','-r','--r',':r'};
leg_hdl = cell(length(M),1);
for i=1:length(M)
    [dt,ebR,ebI,ratio] = returnInformation(rom_data,M(i),t_fom,fom_data.dts);
    leg_hdl{i} = strcat('$M=',strcat(num2str(M(i)),'$'));
    figure(1) %ratio
    plot(t_fom,ratio,linestyle{i},'LineWidth',2);
    hold on
    figure(2) %eb_r
    plot(t_fom,ebR,linestyle{i},'LineWidth',2);
    hold on
    figure(3) %eb_i
    plot(t_fom,ebI,linestyle{i},'LineWidth',2);
    hold on
    figure(4) %phi
    plot(t_fom,atan(ebI./ebR),linestyle{i},'LineWidth',2);
    hold on
end

figure(1)
plot(t_fom(1:end-1),ones(fom_data.n-1,1),'-k','LineWidth',2.5);
hold on  
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Delta t_{\mathrm{ROM},M}/\Delta t_\mathrm{FOM}$','interpreter','latex');
legend(leg_hdl,'interpreter','latex','location','bestoutside');
xlim([0,100]);

figure(2)
plot(t_fom(2:end-1),fom_data.ebR(2:end-1),'-k','LineWidth',2.5);
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\rho(D)$','interpreter','latex');
xlim([0,100]);

figure(3)
plot(t_fom(2:end-1),fom_data.ebI(2:end-1),'-k','LineWidth',2.5);
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\rho(C)$','interpreter','latex');
xlim([0,100]);

figure(4)
plot(t_fom(2:end-1),atan(fom_data.ebI(2:end-1)./fom_data.ebR(2:end-1)),'-k','LineWidth',2.5);
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\varphi$ [rad]','interpreter','latex');
xlim([0,100]);

function [dt,ebR,ebI,ratio] = returnInformation(data,M,t_fom,dt_fom)
rom_data = data.(strcat('rom_',num2str(M)));

t_rom = zeros(rom_data.n,1);
for i=2:rom_data.n
    t_rom(i) = t_rom(i-1)+rom_data.dts(i);
end
dt = interp1(t_rom(1:end-1),rom_data.dts(1:end-1),t_fom,'linear','extrap');
ebR = interp1(t_rom(1:end-1),rom_data.ebR(1:end-1),t_fom,'linear','extrap');
ebI = interp1(t_rom(1:end-1),rom_data.ebI(1:end-1),t_fom,'linear','extrap');
ratio = dt./dt_fom;

end