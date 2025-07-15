function interpol = interpolate(t,vecs,times)

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

interpol = (time2-t)/(time2-time1) * vecs(:,idx_2) ...
         + (t-time1)/(time2-time1) * vecs(:,idx_1);