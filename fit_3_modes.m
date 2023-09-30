function [error_vector] = fit_3_modes(params, alpha, Dmin, Dmax, coff1, coff2, obs1, obs2, obs3, x, sigma1, sigma2, sigma3)
% N0 > 0; mu > -1; lambda > 0
% params = [N0 mu lambda]; input = bins; actual_output = M_sd; x = fit moments

params(2) = max(-1,params(2)); % gammainc requires mu > -1
params(3) = max(0,params(3));
params(1) = max(0,params(1));
params(5) = max(-1,params(5)); % gammainc requires mu > -1
params(6) = max(0,params(6));
params(4) = max(0,params(4));
params(8) = max(-1,params(8)); % gammainc requires mu > -1
params(9) = max(0,params(9));
params(7) = max(0,params(7));

N0 = params(1)/gamma(params(2)+2)*params(3)^2;
mu = params(2);
lambda = params(3);
N02 = params(4)/gamma(params(5)+2)*params(6)^2;
mu2 = params(5);
lambda2 = params(6);
N03 = params(7)/gamma(params(8)+2)*params(9)^2;
mu3 = params(8);
lambda3 = params(9);

fit1 = zeros(size(obs1));
fit2 = zeros(size(obs2));
fit3 = zeros(size(obs3));
f21 = zeros(5,1);
f23 = zeros(5,1);
cfit = zeros(5,1);
loginc1 = zeros(5,1);
loginc2 = zeros(5,1);
loginc3 = zeros(5,1);
loginc4 = zeros(5,1);
loginc5 = zeros(5,1);
loginc6 = zeros(5,1);
loginc7 = zeros(5,1);
loginc8 = zeros(5,1);
loginc9 = zeros(5,1);
loggam1 = zeros(5,1);
loggam2 = zeros(5,1);
loggam3 = zeros(5,1);

for i = 1:5
    loginc1(i) = log( gammainc(lambda*coff1,mu+i) -...
        gammainc(lambda*Dmin,mu+i) );
    loginc2(i) = log( gammainc(lambda*coff2,mu+i) -...
        gammainc(lambda*coff1,mu+i) );
    loginc3(i) = log( gammainc(lambda*Dmax,mu+i) -...
        gammainc(lambda*coff2,mu+i) );
    loginc4(i) = log( gammainc(lambda2*coff1,mu2+i) -...
        gammainc(lambda2*Dmin,mu2+i) );
    loginc5(i) = log( gammainc(lambda2*coff2,mu2+i) -...
        gammainc(lambda2*coff1,mu2+i) );
    loginc6(i) = log( gammainc(lambda2*Dmax,mu2+i) -...
        gammainc(lambda2*coff2,mu2+i) );
    loginc7(i) = log( gammainc(lambda3*coff1,mu3+i) -...
        gammainc(lambda3*Dmin,mu3+i) );
    loginc8(i) = log( gammainc(lambda3*coff2,mu3+i) -...
        gammainc(lambda3*coff1,mu3+i) );
    loginc9(i) = log( gammainc(lambda3*Dmax,mu3+i) -...
        gammainc(lambda3*coff2,mu3+i) );
    loggam1(i) = gammaln(mu+i); % use gammaln to avoid overflow
    loggam2(i) = gammaln(mu2+i);
    loggam3(i) = gammaln(mu3+i);
    f21(i) = N02 ./ ( lambda2.^(i) ) .* exp(loginc4(i)+loggam2(i));
    f23(i) = N02 ./ ( lambda2.^(i) ) .* exp(loginc6(i)+loggam2(i));
    fit1(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc1(i)+loggam1(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loginc4(i)+loggam2(i))...
        + N03 ./ ( lambda3.^(i) ) .* exp(loginc7(i)+loggam3(i));
    fit2(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc2(i)+loggam1(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loginc5(i)+loggam2(i))...
        + N03 ./ ( lambda3.^(i) ) .* exp(loginc8(i)+loggam3(i));
    fit3(i) = N0 ./ ( lambda.^(i) ) .* exp(loginc3(i)+loggam1(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loginc6(i)+loggam2(i))...
        + N03 ./ ( lambda3.^(i) ) .* exp(loginc9(i)+loggam3(i));
    cfit(i) = N0 ./ ( lambda.^(i) ) .* exp(loggam1(i))...
        + N02 ./ ( lambda2.^(i) ) .* exp(loggam2(i))...
        + N03 ./ ( lambda3.^(i) ) .* exp(loggam3(i));
end

loginc10 = loginc1(1);
loginc11 = loginc4(1);
loginc12 = loginc7(1);
loginc13 = loginc2(1);
loginc14 = loginc5(1);
loginc15 = loginc8(1);
loginc16 = loginc3(1);
loginc17 = loginc6(1);
loginc18 = loginc9(1);
loggam4 = loggam1(1); % use gammaln to avoid overflow
loggam5 = loggam2(1);
loggam6 = loggam3(1);
fit4 = N0 ./ ( lambda.^1 ) .* exp(loginc10+loggam4);
fit5 = N02 ./ ( lambda2.^1 ) .* exp(loginc11+loggam5);
fit6 = N03 ./ ( lambda3.^1 ) .* exp(loginc12+loggam6);
fit7 = N0 ./ ( lambda.^1 ) .* exp(loginc13+loggam4);
fit8 = N02 ./ ( lambda2.^1 ) .* exp(loginc14+loggam5);
fit9 = N03 ./ ( lambda3.^1 ) .* exp(loginc15+loggam6);
fit10 = N0 ./ ( lambda.^1 ) .* exp(loginc16+loggam4);
fit11 = N02 ./ ( lambda2.^1 ) .* exp(loginc17+loggam5);
fit12 = N03 ./ ( lambda3.^1 ) .* exp(loginc18+loggam6);
fit13 = N0 .* lambda.^mu .* 0.008.^mu .* exp(-lambda*0.008);
fit14 = N02 .* lambda2.^mu2 .* 0.008.^mu2 .* exp(-lambda2*0.008);

linc1 = log(5*(obs1(1)+obs2(1)+obs3(1))/cfit(1));
linc2 = log(2*(obs1(5)+obs2(5)+obs3(5))/cfit(5));

edt1 = (obs1(2)/obs1(1)-Dmin)/(f21(2)/f21(1)-Dmin)/alpha-1;
edt2 = (f23(4)/f23(3)-Dmin)/(obs3(4)/obs3(3)-Dmin)/alpha-1;
error_vector = zeros(1,25);
for i=1:3
    error_vector(i) = (fit1(1+x(i)) - obs1(1+x(i)))/sigma1(i);
end
for i=1:3
    error_vector(i+3) = (fit2(1+x(i)) - obs2(1+x(i)))/sigma2(i);
end
for i=1:3
    error_vector(i+6) = (fit3(1+x(i)) - obs3(1+x(i)))/sigma3(i);
end
for i=1:3
    error_vector(i+9) = (fit1(1+x(i)) + fit2(1+x(i)) + fit3(1+x(i)) - obs1(1+x(i)) - obs2(1+x(i)) - obs3(1+x(i)))/sqrt(sigma1(i)^2+sigma2(i)^2+sigma3(i)^2);
end

error_vector(13) = 1e3*edt1*(1+sign(edt1));
error_vector(14) = 1e3*edt2*(1+sign(edt2));
error_vector(15) = 1e3*(1+params(2)-params(5))*(1-sign(params(5)-params(2)));
error_vector(16) = 1e3*(1+params(5)-params(8))*(1-sign(params(8)-params(5)));
error_vector(17) = 1e3*(params(6)/params(3))*(1-sign(params(3)-params(6)));
error_vector(18) = 1e3*(params(9)/params(6))*(1-sign(params(6)-params(9)));
error_vector(19) = 1e3*(fit5+fit6)/fit4*(1-sign(alpha^2*fit4/(fit5+fit6)+alpha^2-1));
error_vector(20) = 1e3*(fit7+fit9)/fit8*(1-sign(alpha^2*fit8/(fit7+fit9)+alpha^2-1));
error_vector(21) = 1e3*(fit10+fit11)/fit12*(1-sign(alpha^2*fit12/(fit10+fit11)+alpha^2-1));
error_vector(22) = 1e2*(-linc1)*(1-sign(linc1));
error_vector(23) = 1e2*(-linc2)*(1-sign(linc2));
error_vector(24) = 1e2*(fit14/fit13-alpha^2/(1-alpha^2))*(1-sign(alpha^2*fit13-(1-alpha^2)*fit14));
error_vector(25) = 1e2*((1-alpha^2)*obs1(1)-fit4)/sigma1(1)*(1-sign(fit4-(1-alpha^2)*obs1(1)));
