V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [snapshots.uh_total snapshots.vh_total]';

B = options.rom.B;
basis = B;
Om = options.grid.Om;
F_inhom = options.rom.F_inhom;

%% compute best approximation and error
% V_best = basis*(basis'*(Om.*(snapshots_V_total))); ... % homogeneous part
V_best = basis*(basis'*(Om.*(snapshots_V_total))) ... % homogeneous part
       + F_inhom*((F_inhom'*(Om.*(F_inhom))\(F_inhom'*(Om.*snapshots_V_total)))); % inhomogeneous part  !!!does not work at all!!!

V_ref    = snapshots_V_total;

V_2_ref  = weightedL2norm(V_ref,options.grid.Om);
V_2_ref_avg = sum(V_2_ref)/numel(V_2_ref);

error_V_best = weightedL2norm(V_best - snapshots_V_total,Om)/V_2_ref_avg;



color = 'r';
linestyle = '--';
linewidth = 1;
name = "best approximation error";

figure
plot(t_snapshot,error_V_best,'color',color,'linestyle', linestyle, 'linewidth', linewidth,'displayname',name);
hold on


%% compute ROM error
error_V_ROM = zeros(length(time),1);
for i=1:length(time)
    t = time(i);
    V_FOM = interpolate(t,snapshots_V_total,t_snapshot);
    V_ROM = V_total(:,i);
    error_V_ROM(i) = weightedL2norm(V_FOM-V_ROM,Om)/V_2_ref_avg;
end

linestyle = '-';
name = "ROM error";

plot(time,     error_V_ROM,'color',color,'linestyle', linestyle, 'linewidth', linewidth,'displayname',name);



set(gca,'Yscale','log');
% xlabel('t')
%     set(gcf, 'Position', [100, 100, 400, 300])
% % legend('show','NumColumns',2,'Orientation','horizontal')
legend('show','NumColumns',3,'Orientation','vertical')
%     set(gcf, 'Position', [100, 100, 400, 300])
%     ylim([1e-9 1])
% grid on