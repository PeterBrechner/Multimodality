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
        a = bins-0.5*bins_diff;
        b = bins+0.5*bins_diff;
        M_trap = 1/2 * sum(1e4 * Nd_total .* (a.^(ii)+b.^(ii)) .* bins_diff);
        M(j) = M_trap;
    elseif intmethod == 3
        a = bins-0.5*bins_diff;
        b = bins+0.5*bins_diff;
        M_simpson = 1/6 * sum(1e4 * Nd_total .* (a.^(ii)+4*bins.^(ii)+b.^(ii)) .* bins_diff);
        M(j) = M_simpson;
    elseif intmethod == 4
        a = bins-0.5*bins_diff;
        b = bins+0.5*bins_diff;
        a2b3 = (2*a+b)/3;
        ab23 = (a+2*b)/3;
        M_simps38 = 1/8 * sum(1e4 * Nd_total .* (a.^(ii)+3*a2b3.^(ii)+3*ab23.^(ii)+b.^(ii)) .* bins_diff);
        M(j) = M_simpson38;
    elseif intmethod == 5
        a = Nd_total(1:end-1);
        b = Nd_total(2:end);
        c = 1/6*bins_diff(1:end-1);
        d = 1/6*bins_diff(2:end);
        e = bins(2:end)-0.25*bins_diff(2:end)-0.25*bins_diff(1:end-1);
        f = bins(1)-0.5*bins_diff(1);
        g = bins(end)+0.5*bins_diff(end);
        mult_vector = ones(1,length(bins)-1);
        mult_vector(1:2:length(bins)-1) = 4*mult_vector(1:2:length(bins)-1);
        mult_vector(2:2:length(bins)-1) = 2*mult_vector(2:2:length(bins)-1);
        M_low = 1e4 * bins_diff(1)/3 * Nd_total(1) * f^ii;
        M_high = 1e4 * bins_diff(end)/3 * Nd_total(end) * g^ii;
        M_simpscomposite = M_low + M_high + 1e4 * sum(mult_vector.*(c.*a+d.*b).*e.^ii);
        M(j) = M_simpscomposite;
    end
end
end
end
