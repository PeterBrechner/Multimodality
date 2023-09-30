function [error_vector] = fit_bimodal2(params, alpha, Dmin, Dmax, coff, obs1, obs2, x, sigma1, sigma2)
% N0 > 0; mu > -1; lambda > 0
% params = [N0 mu lambda]; input = bins; actual_output = M_sd; x = fit moments

params(1) = max(0,params(1));
params(2) = max(-1,params(2));
params(3) = max(0,params(3));
params(4) = max(0,params(4));
params(5) = max(-1,params(5));
params(6) = max(0,params(6));

N0 = params(1)/gamma(params(2)+2)*params(3)^2;
mu = params(2);
lambda = params(3);
N02 = params(4)/gamma(params(5)+2)*params(6)^2;
mu2 = params(5);
lambda2 = params(6);

fit1 = zeros(size(obs1));
fit2 = zeros(size(obs2));
f21 = zeros(5,1);
f23 = zeros(5,1);
cfit = zeros(5,1);
loginc1 = zeros(5,1);
loginc2 = zeros(5,1);
loginc3 = zeros(5,1);
loginc4 = zeros(5,1);
loggam = zeros(5,1);
loggam2 = zeros(5,1);

for i = 1:5
    loginc1(i) = log( gammainc(lambda*coff,mu+i) -...
        gammainc(lambda*Dmin,mu+i) );
    loginc2(i) = log( gammainc(lambda*Dmax,mu+i) -...
        gammainc(lambda*coff,mu+i) );
    loginc3(i) = log( gammainc(lambda2*coff,mu2+i) -...
        gammainc(lambda2*Dmin,mu2+i) );
    loginc4(i) = log( gammainc(lambda2*Dmax,mu2+i) -...
        gammainc(lambda2*coff,mu2+i) );
    loggam(i) = gammaln(mu+i); % use gammaln to avoid overflow
    loggam2(i) = gammaln(mu2+i);
    f21(i) =  N02 ./ ( lambda2.^(i) ) .* exp(loginc3(i)+loggam2(i));
    f23(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc2(i)+loggam(i));
    fit1(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc1(i)+loggam(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loginc3(i)+loggam2(i));
    fit2(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc2(i)+loggam(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loginc4(i)+loggam2(i));
    cfit(i) = N0 ./ ( lambda.^(i) ) .* exp(loggam(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loggam2(i));
end
loginc5 = loginc1(1);
loginc6 = loginc3(1);
loginc7 = loginc2(1);
loginc8 = loginc4(1);
loggam3 = loggam(1); % use gammaln to avoid overflow
loggam4 = loggam2(1);
fit3 = N0 ./ ( lambda.^1 ) .* exp(loginc5+loggam3);
fit4 = N02 ./ ( lambda2.^1 ) .* exp(loginc6+loggam4);
fit5 = N0 ./ ( lambda.^1 ) .* exp(loginc7+loggam3);
fit6 = N02 ./ ( lambda2.^1 ) .* exp(loginc8+loggam4);
fit7 = N0 .* lambda.^mu .* 0.03.^mu .* exp(-lambda*0.03);
fit8 = N02 .* lambda2.^mu2 .* 0.03.^mu2 .* exp(-lambda2*0.03);

linc1 = log(5*(obs1(1)+obs2(1))/cfit(1));
linc2 = log(2*(obs1(5)+obs2(5))/cfit(5));

edt1 = (obs1(2)/obs1(1)-Dmin)/(f21(2)/f21(1)-Dmin)/alpha-1;
edt2 = (f23(4)/f23(3)-Dmin)/(obs2(4)/obs2(3)-Dmin)/alpha-1;

error_vector = zeros(1,19);
for i=1:3
    error_vector(i) = (fit1(1+x(i)) - obs1(1+x(i)))/sigma1(i);
end
for i=1:3
    error_vector(i+3) = (fit2(1+x(i)) - obs2(1+x(i)))/sigma2(i);
end
for i=1:3
    error_vector(i+6) = (fit1(1+x(i)) + fit2(1+x(i)) - obs1(1+x(i)) - obs2(1+x(i)))/sqrt(sigma1(i)^2+sigma2(i)^2);
end
error_vector(10) = 1000*edt1*(1+sign(edt1));
error_vector(11) = 1000*edt2*(1+sign(edt2));
error_vector(12) = 1000*(1+params(2)-params(5))*(1-sign(params(5)-params(2)));
error_vector(13) = 1000*(params(6)/params(3))*(1-sign(params(3)-params(6)));
error_vector(14) = 1000*fit4/fit3*(1-sign(alpha^2*fit3/fit4+alpha^2-1));
error_vector(15) = 1000*fit5/fit6*(1-sign(alpha^2*fit6/fit5+alpha^2-1));
error_vector(16) = 100*(-linc1)*(1-sign(linc1));
error_vector(17) = 100*(-linc2)*(1-sign(linc2));
error_vector(18) = 100*(fit8/fit7-alpha^2/(1-alpha^2))*(1-sign(alpha^2*fit7-(1-alpha^2)*fit8));
error_vector(19) = 100*((1-alpha^2)*obs1(1)-fit3)/sigma1(1)*(1-sign(fit3-(1-alpha^2)*obs1(1)));
