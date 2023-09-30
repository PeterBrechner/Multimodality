function [tr, ba, bb, un, hix, hix2, mp3, mp4, mp5, mp6, upr] =...
    TrimodalityTest(uonly, testsH, testsL, alpha, time, intmethod, fit_moments,...
    bins, bins_diff, sd, iwc2, rawcount, error, decider, conf, cts, lamlam, hdb2)


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
%hdb2 = control for whether habits are used to set testsH

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
if hdb2
    b = testsH;
end

%Convert to old definition of alpha for coding purposes
alpha = sqrt(1-alpha);

%Calculate maximum D among non-outlier counts for each SD
upr = zeros(sz(1),1);
for j=1:sz(1)
    upr(j) = calcMax(sd(j,:),bins,bins_diff);
end

%Run tests

if uonly
    un = 1:sz(1);
    mp6 = Fit1Mode(vecMin(un), vecMax(un), intmethod, sd(un,:),...
        bins, bins_diff, fit_moments, error(un,:), upr(un));
    return
end

if ~hdb2
    %Find reference fit for first large diameter dip test
    disp(strcat("Starting Testing for Large Mode Dcutoff >= 1 mm: ", string(datetime("now"))));
    [minparams2] = FitCenterMode(vecMin,vecMax,intmethod,sd,bins,bins_diff,0,error,upr);    
    %Estimate large diameter cutoff assuming no small diameter bimodality
    %Test only larger cutoffs to avoid misclassifying bimodal1 as bimodal2
    [~, minr2, b, minparams] = findCutoff(alpha, testsH(4:end), 0.03, b, minparams,...
        minparams2, 0, intmethod, sd, bins, bins_diff, error, rawcount,...
        decider, conf, cts, lamlam, max(testsL), upr);
    %minparamsx = minparams;
    hi2 = b;
end

%Find reference fit for small diameter dip test
disp(strcat("Starting Testing for Small Mode: ", string(datetime("now"))));
[minparams2] = FitCenterMode(vecMin,b,intmethod,sd,bins,bins_diff,1,error,upr);
%Estimate small diameter cutoff
[~, minr, b, minparams] = findCutoff(alpha, testsL, 0.03, b, minparams,...
    minparams2, 1, intmethod, sd, bins, bins_diff, error, rawcount,...
    decider, conf, cts, lamlam, min(testsH), upr);
hi = b;


%Find reference fit for final large diameter dip test
disp(strcat("Starting Testing for Large Mode: ", string(datetime("now"))));
[minparams2] = FitCenterMode(b,vecMax,intmethod,sd,bins,bins_diff,0,error,upr);
if ~hdb2
    %minparams = minparamsx;
    %[coff2, minr2, b, minparams] = findCutoff(alpha, testsH, 0.03, b, minparams,...
    %    minparams2, 0, intmethod, sd, bins, bins_diff, error, rawcount,...
    %    decider, conf, cts, lamlam, max(testsL), upr);
    %minparamsx = minparams;
    %hi2 = b
    
    %Estimate large diameter cutoff when small mode exists
    smex = find((b > 0) | (hi2 ~= sz(2)));
    nsmex = find((b == 0) & (hi2 == sz(2)));
    %Consider smaller large diameter cutoffs now that bimodal1 is tested
    %minparams(nsmex,:) = minparamsx(nsmex,:);
    disp(strcat("Correcting Large Mode Dcutoff for Currently Multimodal PSDs: ", string(datetime("now"))));
    [~, minr2(smex), b(smex), minparams(smex,:)] = findCutoff(alpha, ...
        testsH, 0.03, b(smex), minparams(smex,:), minparams2(smex,:), ...
        0, intmethod, sd(smex,:), bins, bins_diff, error(smex,:), ...
        rawcount(smex,:), decider, conf, cts, lamlam, max(testsL), upr(smex));
    disp(strcat("Starting Testing for Large Mode Dcutoff < 1 mm for Currently Unimodal PSDs: ", string(datetime("now"))));
    [~, minr2(nsmex), b(nsmex), minparams(nsmex,:)] = findCutoff(alpha, ...
        testsH(1:3), 0.03, b(nsmex), minparams(nsmex,:), minparams2(nsmex,:), ...
        0, intmethod, sd(nsmex,:), bins, bins_diff, error(nsmex,:), ...
        rawcount(nsmex,:), decider, conf, cts, lamlam, max(testsL), upr(nsmex));
    %hi2(smex) = b(smex);
    hi2 = b;
    
else
    minr2 = zeros(sz(1),1);
    coff2 = zeros(sz(1),1);
    for j=1:sz(1)
        testsHxx = bins(testsH(j))+0.5*bins_diff(testsH(j));
        %testsHx = [bins(testsH(j)-1)+0.5*bins_diff(testsH(j)-1),...
        %    bins(testsH(j))+0.5*bins_diff(testsH(j)),...
        %    bins(testsH(j)+1)+0.5*bins_diff(testsH(j)+1)]
        testsHx = [0.85001*testsHxx, testsHxx, 1.14999*testsHxx];
        [coff2(j), minr2(j), b(j), minparams(j,:)] = findCutoff(alpha, ...
            testsHx, 0.03, b(j), minparams(j,:), minparams2(j,:), 0, ...
            intmethod, sd(j,:), bins, bins_diff, error(j,:), ...
            rawcount(j,:), decider, conf, cts, lamlam, max(testsL), upr(j));
    end
    hi2 = b;
end

%Fill cutoff bin number arrays
%hi2 = b;
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
%mp3 = zeros(length(tr),9);
%mp4 = zeros(length(ba),6);
%mp5 = zeros(length(bb),6);
%mp6 = zeros(length(un),3);

%Fit multimodal distributions
disp(strcat("Fitting Trimodal PSDs: ", string(datetime("now"))));
mp3 = Fit3Modes(alpha, mag(tr,:), hix(tr), hix2(tr), intmethod, sd(tr,:),...
    bins, bins_diff, fit_moments, error(tr,:), upr(tr));
disp(strcat("Fitting Bimodal1 PSDs: ", string(datetime("now"))));
mp4 = FitAnyModes(alpha, mag(ba,:), hix(ba), hix2(ba), intmethod, sd(ba,:),...
    bins, bins_diff, fit_moments, error(ba,:), upr(ba));
disp(strcat("Fitting Bimodal2 PSDs: ", string(datetime("now"))));
mp5 = FitAnyModes(alpha, mag(bb,:), hix(bb), hix2(bb), intmethod, sd(bb,:),...
    bins, bins_diff, fit_moments, error(bb,:), upr(bb));
disp(strcat("Fitting Unimodal PSDs: ", string(datetime("now"))));
mp6 = Fit1Mode(vecMin(un), vecMax(un), intmethod, sd(un,:),...
    bins, bins_diff, fit_moments, error(un,:), upr(un));
%PlotMedians2(1e4*bins, 1e3*sd(tr,:), 1e3*sd(ba,:), 1e3*sd(bb,:), 1e3*sd(un,:), tr, ba, bb,...
% un, iwc2, "g/m^3", 'HIWCimagesLatest/guess');

end
