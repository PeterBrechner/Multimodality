function [error_vector] = fit_1_mode(params, Dmin, Dmax, coff1, coff2, obs, x, sigma)
% N0 > 0; mu > -1; lambda > 0
% params = [N0 mu lambda]; input = bins; actual_output = M_sd; x = fit moments

params(1) = max(0,params(1));
params(2) = max(-1,params(2));
params(3) = max(0,params(3));

N0 = params(1)/gamma(params(2)+2)*params(3)^2;
mu = params(2);
lambda = params(3);

fit = zeros(size(obs));
cfit = zeros(3,1);
loginc = zeros(3,1);

for i = 1:length(x)
    loginc(i) = log( gammainc(lambda*coff2,1+mu+x(i)) -...
        gammainc(lambda*coff1,1+mu+x(i)) );
    loggam = gammaln(1+mu+x(i)); % use gammaln to avoid overflow
    fit(i) = N0 ./ ( lambda.^(1+x(i)) ) .* exp(loginc(i)+loggam);
    cfit(i) = N0 ./ ( lambda.^(1+x(i)) ) .* exp(loggam);
end

linc1 = log(5*obs(1)/cfit(1));
linc2 = log(2*obs(3)/cfit(3));

error_vector = zeros(1,5);
for i=1:3
    error_vector(i) = (fit(i) - obs(i))/sigma(i);
end
error_vector(4) = 100*(-linc1)*(1-sign(linc1));
error_vector(5) = 100*(-linc2)*(1-sign(linc2));
