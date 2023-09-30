function overlap = overlapTest(Dmin, Dmax, N0, mu, lambda, alpha, sd, bins,...
    bins_diff, sqError, conf, ii, ii2)

%Output: overlap measures the maximum fraction of the iith moment...
%    over edge bins that can be explained by the center mode

%idx1 denotes the edge diameter bins
%    Dmin is the large diameter cutoff
%    Or Dmax is the small diameter cutoff
idx1 = find(bins > Dmin & bins < Dmax);

%Calculate the maximum fraction of the iith moment over edge bins...
%    that can be explained by the center mode
loginc = log( gammainc(lambda*Dmax,1+mu+ii) -...
    gammainc(lambda*Dmin,1+mu+ii) );
loggam = gammaln(1+mu+ii); % use gammaln to avoid overflow
%tol = uncertainty
tol = sqrt(chi2inv(conf,3))*calcError(sqError(idx1), bins(idx1), ii);
%For bimodality, need obs-tol > fit/alpha
%    Define fit=alpha*(tol+fit/alpha) for mathematical simplicity
fit = N0 ./ ( lambda.^(1+mu+ii) ) .* exp(loginc+loggam) + alpha*tol;
%Calculate the observed iith moment
obs = intMethods(1, sd(idx1), bins(idx1), bins_diff(idx1), ii);

minD = bins(1)-0.5*bins_diff(1);
obs2 = intMethods(1, sd(idx1), bins(idx1), bins_diff(idx1), ii2+1);
obs1 = intMethods(1, sd(idx1), bins(idx1), bins_diff(idx1), ii2);
obs3 = obs2/obs1-minD;
loginc1 = log( gammainc(lambda*Dmax,1+mu+ii2) -...
    gammainc(lambda*Dmin,1+mu+ii2) );
loggam1 = gammaln(1+mu+ii2);
loginc2 = log( gammainc(lambda*Dmax,2+mu+ii2) -...
    gammainc(lambda*Dmin,2+mu+ii2) );
loggam2 = gammaln(2+mu+ii2);
fit3 = N0 ./ ( lambda.^(2+mu+ii2) ) .* exp(loginc2+loggam2) ./...
    (N0 ./ ( lambda.^(1+mu+ii2) ) .* exp(loginc1+loggam1));
fit3 = fit3-minD;
shoulderAlpha = fit3./obs3;
if minD == Dmin
    shoulderAlpha = 1/shoulderAlpha;
end

overlap = max(sqrt(fit./obs), shoulderAlpha);

end

