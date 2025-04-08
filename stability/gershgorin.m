function bound = gershgorin(M)
  absM = abs(M);
  bound = max(sum(absM,2));
end
