function [tr, ba, bb, un, hix, hix2, mp3, mp4, mp5, mp6, upr] =...
    TrimodalityTest(uonly, testsH, testsL, alpha, time, intmethod, fit_moments,...
    bins, bins_diff, sd, iwc2, rawcount, error, decider, conf, cts, lamlam);


%Tests whether distributions are unimodal, bimodal, or trimodal

%Outputs:
%tr = indices of trimodal distributions
%ba = indices of distributions with small diameter bimodality
%bb = indices of distributions with large diameter bimodality
%un = indices of unimodal distributions
%hix = highest bin number below lower cutoff
%hix2 = highest bin number below upper cutoff
%mp3 = fit parameters for trimodal distributions
%mp4 = fit parameters for distributions with small diameter bimodality
%mp5 = fit parameters for distributions with large diameter bimodality
%mp6 = fit parameters for unimodal distributions

%Inputs:
%uonly = flag to force unimodal distributions
%testsH = cutoffs to test for large diameter bimodality, in ascending order
%testsL = cutoffs to test for small diameter bimodality, in descending order
%alpha = "sensitivity"
%time = time of size distribution (s)
%intmethod = integration method used by TrimodalityTest.m
%fit_moments = fitting moments used
%bins = bins of size distribution (cm)
%bins_diff = bin widths (cm)
%sd = size distribution (cm^-3 um^-1)
%iwc2 = ice water content (g m^-3)
%rawcount = raw counts in each bin
%error = squared error in bin_concentration (cm^-6)
%decider = method of choosing cutoffs
%conf = confidence level
%cts = minimum counts in center mode
%lamlam = estimated ratio in lambda between consecutive modes

%Initialize outputs
sz = size(sd);
tr = time;
ba = time;
bb = time;
un = time;
hix = 0;
hix2 = 0;
mp3 = 0;
mp4 = 0;
mp5 = 0;
mp6 = 0;
if sz(1) == 0
    return
end

%Initialize arrays
hix = zeros(sz(1),1);
hix2 = sz(2)*ones(sz(1),1);
minparams = zeros(sz(1),3);
vecMin = hix;
vecMax = hix2;
b = zeros(sz(1),1);

%Convert to old definition of alpha for coding purposes
alpha = sqrt(1-alpha);

%Calculate maximum D among non-outlier counts for each SD
upr = zeros(sz(1),1);
for j=1:sz(1)
    upr(j) = calcMax(sd(j,:),bins,bins_diff,alpha);
end

%Run tests

if uonly
    un = 1:sz(1);
    mp6 = Fit1Mode(alpha, vecMin(un), vecMax(un), intmethod, sd(un,:),...
        bins, bins_diff, fit_moments, error(un,:), upr(un));
    return
end

%Find reference fit for first large diameter dip test
[minparams2] = FitCenterMode(alpha,vecMin,vecMax,intmethod,sd,bins,bins_diff,0,error,upr);

%Estimate large diameter cutoff assuming no small diameter bimodality
[coff2, minr2, b, minparams] = findCutoff(alpha, testsH, 0.03, b, minparams,...
    minparams2, 0, intmethod, sd, bins, bins_diff, error, iwc2, rawcount,...
    decider, conf, cts, lamlam, max(testsL), upr);

%Find reference fit for small diameter dip test
hi2 = b
[minparams2] = FitCenterMode(alpha,vecMin,b,intmethod,sd,bins,bins_diff,1,error, upr);

%Estimate small diameter cutoff
[coff, minr, b, minparams] = findCutoff(alpha, testsL, 0.03, b, minparams,...
    minparams2, 1, intmethod, sd, bins, bins_diff, error, iwc2, rawcount,...
    decider, conf, cts, lamlam, min(testsH), upr);

%Find reference fit for final large diameter dip test
hi = b
[minparams2] = FitCenterMode(alpha,b,vecMax,intmethod,sd,bins,bins_diff,0,error, upr);

%Estimate large diameter cutoff
[coff2, minr2, b, minparams] = findCutoff(alpha, testsH, 0.03, b, minparams,...
    minparams2, 0, intmethod, sd, bins, bins_diff, error, iwc2, rawcount,...
    decider, conf, cts, lamlam, max(testsL), upr);

%Fill cutoff bin number arrays
hi2 = b
hix = hi;
hix2 = hi2;

%Fill modality index arrays
tr = find(minr < 0 & minr2 < 0);
ba = find(minr < 0 & minr2 >= 0);
bb = find(minr >= 0 & minr2 < 0);
un = find(minr >= 0 & minr2 >= 0);

%Calculate fit parameters used to initialize multimodal fits
mpx = minparams(:,2:3);
mag = intMethods(intmethod, sd, bins, bins_diff, 1);
mag = [mag, mpx];

%Initialize fit parameter arrays
mp3 = zeros(length(tr),9);
mp4 = zeros(length(ba),6);
mp5 = zeros(length(bb),6);
mp6 = zeros(length(un),3);

%Fit multimodal distributions
mp3 = Fit3Modes(alpha, mag(tr,:), hix(tr), hix2(tr), intmethod, sd(tr,:),...
    bins, bins_diff, fit_moments, error(tr,:), upr(tr));
mp4 = FitAnyModes(alpha, mag(ba,:), hix(ba), hix2(ba), intmethod, sd(ba,:),...
    bins, bins_diff, fit_moments, error(ba,:), upr(ba));
mp5 = FitAnyModes(alpha, mag(bb,:), hix(bb), hix2(bb), intmethod, sd(bb,:),...
    bins, bins_diff, fit_moments, error(bb,:), upr(bb));
mp6 = Fit1Mode(alpha, vecMin(un), vecMax(un), intmethod, sd(un,:),...
    bins, bins_diff, fit_moments, error(un,:), upr(un));
PlotMedians2(1e4*bins, 1e3*sd(tr,:), 1e3*sd(ba,:), 1e3*sd(bb,:), 1e3*sd(un,:), tr, ba, bb,...
 un, iwc2, "g/m^3", 'HIWCimagesLatest/guess');

end
