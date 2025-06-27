%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy
ken_hdl.kdecay = (k(1:n-1)-k(1))/k(1);
ken_hdl.t = time(1:n-1);

figure(1)
plot(ken_hdl.t,ken_hdl.kdecay)
grid on
hold on
xlabel('$t$','interpreter','latex');
ylabel('$(k(t)-k(0))/k(0)$','interpreter','latex');


%% momentum
figure
% plot(time,umom,'s-')
% hold on
% plot(time,vmom,'o-')
if (options.rom.rom == 1 )
    umom0 = snapshots.umom(1);
else
    umom0 = umom(1);
end
% open('results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');
hold on
semilogy(time,abs(umom-umom0)/umom0,'s--')
title('u-momentum error');
% saveas(gcf,'results/shear_layer_ROM/mom_error_ROM_inviscid_momcons.fig');

%% vorticity
% compare e.g. with PhD thesis figure 3.4

figure
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx+1,Npy+1);
% for Re=1000: labels = -4:0.5:4;
labels= 20;
contour(x,y,omega',labels);
axis square
% colorbar
grid