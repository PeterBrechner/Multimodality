function [M1] = find_starting_M1(mu1, lc1, mu2, lc2)

M1 = gamma(2+mu1)/gamma(2+mu2)*(lc2)^(2+mu2)/(lc1)^(2+mu1)*exp(lc1-lc2);

end
