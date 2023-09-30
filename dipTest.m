function [ratio, overlap] = dipTest(idx1, idx2, idx3, ll, ml, ul, uu, N01,...
    mu1, lambda1, conf, alpha, overlap, sd, bins, bins_diff, sqError, ii)

%Outputs:
%ratio is the negative of the depth of the dip
%    beyond the threshold for bimodality
%overlap is the output of the overlap test for the SD being tested

%Calculate observed partial moments
lower = intMethods(1, sd(idx1), bins(idx1), bins_diff(idx1), ii);
mid = intMethods(1, sd(idx2), bins(idx2), bins_diff(idx2), ii);
upper = intMethods(1, sd(idx3), bins(idx3), bins_diff(idx3), ii);

%Calculate normalized fit partial moments
L = gammainc(ml*lambda1,1+ii+mu1)-gammainc(ll*lambda1,1+ii+mu1);
M = gammainc(ul*lambda1,1+ii+mu1)-gammainc(ml*lambda1,1+ii+mu1);
U = gammainc(uu*lambda1,1+ii+mu1)-gammainc(ul*lambda1,1+ii+mu1);

%Calculate scaling factor for normalized fit partial moments
mag = N01*gamma(1+ii+mu1)/lambda1^(1+ii+mu1);

%Calculate log(observed partial moment ratio)
ratio = log(mid/sqrt(lower*upper));

%Calculate standard deviations of logs of observed partial moments
errorL = calcError(sqError(idx1), bins(idx1), ii)/lower;
%Modify standard deviation calculation for log of middle partial moment
%    to ensure mid*exp(sqrt(chi2inv(conf,3))*errorM) decreases
%    when mid decreases
errorM = log(1+sqrt(chi2inv(conf,3))...
               *calcError(sqError(idx2), bins(idx2), ii)/mid)...
         /sqrt(chi2inv(conf,3));
%If middle partial moment is 0, modes are perfectly separated
%Set errorM to 0 in this case so SD can pass dip test
if mid == 0
    errorM = 0;
end
errorU = calcError(sqError(idx3), bins(idx3), ii)/upper;

%Calculate tolerances for fit partial moment ratio
%    Tolerances add equivalent of about 1/2 count to each fit partial moment
%    to ensure fit partial moment ratio does not blow up at small and large D
tolL = calcError(sqError(idx1), bins(idx1), ii)/mag;
tolM = calcError(sqError(idx2), bins(idx2), ii)/mag;
tolU = calcError(sqError(idx3), bins(idx3), ii)/mag;

%Calculate log(alpha times fit partial moment ratio with tolerances)
mult = log(alpha*sqrt(M^2+tolM^2)/...
    sqrt(sqrt(L^2+tolL^2)*sqrt(U^2+tolU^2)));

%If middle partial moment is 0, ensure SD passes dip test
%    if other partial moments have enough counts 
if mid == 0
    ratio = mult-1;
end

%Calculate negative depth of dip beyond threshold for bimodality
ratio = ratio-mult+sqrt(chi2inv(conf,3))*...
    sqrt(errorM^2+0.25*errorL^2+0.25*errorU^2);

%If dip test fails, force overlap test to fail
if ratio >= 0
    overlap = 1+overlap;
end
