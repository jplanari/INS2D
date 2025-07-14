function X_bc = get_X_bc(ts,options)

X_bc = zeros(length(get_bc_vector_yBC(0,options)),length(ts));
for jj=1:length(ts)
    t_j = ts(jj);
    X_bc(:,jj) = get_bc_vector_yBC(t_j,options);
end