%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy

ken_hdl.kdecay = (k(1:n)-k(1))/k(1);
ken_hdl.t = time(1:n);

figure(1)
plot(ken_hdl.t,ken_hdl.kdecay)
grid on
hold on
xlabel('$t$','interpreter','latex');
ylabel('$(k(t)-k(0))/k(0)$','interpreter','latex');

%% vorticity
% compare e.g. with PhD thesis figure 3.4

figure
omega = get_vorticity(V,t,options);
omega = reshape(omega,Npx+1,Npy+1);
% for Re=1000: labels = -4:0.5:4;
labels= 20;
contour(x,y,omega',labels);
axis square
cb_handle = colorbar;
grid on
hold on

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
ylabel(cb_handle, '$\omega$','interpreter','latex') 

drawnow % Ensure positions are available

% Shrink the axes slightly to make room for the colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1)+0.05;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.06; % leave space for colorbar
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left, bottom, ax_width, ax_height];
