function bound = gershgorin(A)
      absM = abs(A);
      bound = max(sum(absM,2));
end
