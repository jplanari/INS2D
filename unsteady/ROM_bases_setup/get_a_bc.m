function a_bc = get_a_bc(t,options)

times = options.rom.a_BCs_times;
a_BCs = options.rom.a_BCs;

idx_1 = find(times > t, 1, 'first');
idx_2 = find(times < t, 1, 'last');

if (isempty(idx_2))
    idx_2 = 1;
end

if (isempty(idx_1))
    idx_1 = length(t_snp);
end

time1 = times(idx_1);
time2 = times(idx_2);

a_bc = (time2-t)/(time2-time1) * a_BCs(:,idx_2) ...
     + (t-time1)/(time2-time1) * a_BCs(:,idx_1);