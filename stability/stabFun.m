function R = stabFun(s,z)
    R=1;
    for i=1:s
        R = R + 1/factorial(i)*z^i;
    end
end
