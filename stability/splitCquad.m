function [Cr_tensor] = splitCquad(Cquad,M)
    Cr_tensor = zeros(M,M,M);
    for i=1:M
        cols = i:M:M^2-M+i;
        Cr_tensor(:,:,i) = Cquad(:,cols);
    end
end

