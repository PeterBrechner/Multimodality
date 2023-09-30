function M = intMethods(intmethod, sd, bins, bins_diff, ii)

sd_size = size(sd);
M = zeros(sd_size(1),1);
if intmethod == 0
    musize = size(ii);
    M = zeros(musize);
    Nd_total = sd;
    Nd_total(isnan(Nd_total) == 1) = 0;
    Nd_total(isinf(Nd_total) == 1) = 0;
    for i = 1:length(Nd_total);
        M = M + 1e4 * Nd_total(i) * bins(i).^(ii) * bins_diff(i);
    end
else
for j = 1:sd_size(1)

    Nd_total = sd(j,:); % Nd_total: current sd
    Nd_total(isnan(Nd_total) == 1) = 0;
    Nd_total(isinf(Nd_total) == 1) = 0;
%    [num2str(j),' of ',num2str(sd_size(1))]

    if intmethod == 1
        M(j) = sum(1e4 * Nd_total .* bins.^(ii) .* bins_diff); % M_rect
    elseif intmethod == 2
        M_trap = 1e4 * trapz(bins, Nd_total .* bins.^(ii));
        M(j) = M_trap;
    elseif intmethod == 3
        j1 = 1:length(bins)-1;
        j2 = 2:length(bins);

        a = bins(j1);
        b = bins(j2);
        fa = Nd_total(j1);
        fb = Nd_total(j2);
        ab2 = (a+b)/2;
        yi = interp1(bins, Nd_total, ab2,'pchip');
        simps = (b-a)/6 .* (fa + 4.*yi + fb);

        M_simpson = sum(1e4*simps.*ab2.^ii);
        M(j) = M_simpson;
    elseif intmethod == 4
        j1 = 1:length(bins)-1;
        j2 = 2:length(bins);

        a = bins(j1);
        b = bins(j2);
        fa = Nd_total(j1);
        fb = Nd_total(j2);
        ab2 = (a+b)/2;

        a2b3 = (2*a+b)/3;
        ab23 = (a+2*b)/3;
        combined = sort([a2b3 ab23]);
        yia2b3 = interp1(bins,Nd_total,combined,'pchip');

        yi1 = yia2b3(1:2:end);
        yi2 = yia2b3(2:2:end);

        simps38 = (b-a)/8 .* (fa' + 3.*yi1' + 3 .* yi2' + fb');

        M_simpson38 = 1e4 * sum(simps38.*ab2.^ii);
        M(j) = M_simpson38;
    elseif intmethod == 5
        mult_vector = ones(1,length(bins));
        mult_vector(2:2:length(bins)-1) = 4*mult_vector(2:2:length(bins)-1);
        mult_vector(3:2:length(bins)-1) = 2*mult_vector(3:2:length(bins)-1);
        M_simpscomposite = 1e4 * sum(bins_diff/3.*mult_vector.*Nd_total .*bins.^ii);
        M(j) = M_simpscomposite;
    end
end
end
end
