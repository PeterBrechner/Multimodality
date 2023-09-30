function [sd, sqError, modality, N0s, mus, lambdas, meds] = biasTest(moments, ratios, mus, lambdas, alpha, binEdges, total)

%tests for biases in modality method

%moments = moments where neighboring modes have equal magnitude
	%example: [0, 4] means moment 0 is the same for the small and center modes
	%and moment 4 is the same for the center and large modes
%mus = set of theoretical values for mu
%lambdas = set of theoretical values for lambda
%alpha = sensitivity parameter used for modality
%total = total number of counts allowed
%sd = size distribution
%sqError = squared error in each bin
%modality = theoretical modality type
%N0s = theoretical N0 set

mus = sort(mus);
lambdas = flip(sort(lambdas));
meds = gammaincinv(0.5,1+mus)./lambdas;
total1 = ratios(1)*(lambdas(1)/lambdas(2))^moments(1)*gamma(1+mus(1))/gamma(1+mus(2))...
        *gamma(1+moments(1)+mus(2))/gamma(1+moments(1)+mus(1));
moment41 = (lambdas(1)/lambdas(2))^4*gamma(1+mus(1))/gamma(1+mus(2))...
        *gamma(1+4+mus(2))/gamma(1+4+mus(1));
total3 = ratios(2)*(lambdas(3)/lambdas(2))^moments(2)*gamma(1+mus(3))/gamma(1+mus(2))...
        *gamma(1+moments(2)+mus(2))/gamma(1+moments(2)+mus(3));
moment43 = (lambdas(3)/lambdas(2))^4*gamma(1+mus(3))/gamma(1+mus(2))...
        *gamma(1+4+mus(2))/gamma(1+4+mus(3));
denominator = 1+total1+total3;
total1 = round(total1*total/denominator)
total2 = round(total/denominator)
total3 = round(total3*total/denominator)
if meds(1) <= meds(2)*alpha/(2-alpha) & meds(2) <= meds(3)*alpha/(2-alpha)...
    & total1 >= total2*(1-alpha)*meds(1)/meds(2)...
    & total1 <= total2*moment41/(1-alpha)...
    & total3 >= total2*moment43*(1-alpha)...
    & total3 <= total2*meds(3)/meds(2)/(1-alpha)...
    & total1 >= chi2inv(0.95,3) & total2 >= chi2inv(0.95,3)...
    & total3 >= chi2inv(0.95,3)
	modality = "trimodal";
elseif meds(1) <= meds(2)*alpha/(2-alpha)...
    & total1 >= total2*(1-alpha)*meds(1)/meds(2)...
    & total1 <= total2*moment41/(1-alpha)...
    & (total2 >= chi2inv(0.95,3) | total3 >= chi2inv(0.95,3))...
    & total1 >= chi2inv(0.95,3)
	modality = "bimodal1";
elseif meds(2) <= meds(3)*alpha/(2-alpha)...
    & total3 >= total2*moment43*(1-alpha)...
    & total3 <= total2*meds(3)/meds(2)/(1-alpha)...
    & total2 >= chi2inv(0.95,3) & total3 >= chi2inv(0.95,3)
	modality = "bimodal2";
else
	modality = "unimodal";
end
N0s = [total1*lambdas(1)^(1+mus(1))/gamma(1+mus(1)),...
	total2*lambdas(2)^(1+mus(2))/gamma(1+mus(2)),...
	total3*lambdas(3)^(1+mus(3))/gamma(1+mus(3))];
d = gamrnd(1+mus(1), 1/lambdas(1), total1, 1);
d = [d; gamrnd(1+mus(2), 1/lambdas(2), total2, 1)];
d = [d; gamrnd(1+mus(3), 1/lambdas(3), total3, 1)];
d = d(find(d > binEdges(1) & d < binEdges(end)));
h = histogram(d, binEdges);
sd = h.Values;
binUpper = binEdges(2:end);
binLower = binEdges(1:end-1);
binWidth = binUpper-binLower;
sd = sd./(1e4*binWidth);
sqError = h.Values;
