function [tr, ba, bb, un, hix, hix2, mp3, mp4, mp5, mp6] =...
    Input(time, fit_moments, bins, bins_diff, sd, iwc, rawcount, sqError, idx)


%Call this function to run TrimodalityTest.m


%Outputs of Input.m

%tr = indices of trimodal distributions
    %for example, if 4 is in tr, the fourth size distribution is trimodal
%ba = indices of distributions with small diameter bimodality
%bb = indices of distributions with large diameter bimodality
%un = indices of unimodal distributions
%hix = highest bin number below lower cutoff
    %for example, if hix = 10, there are 10 bins below the lower cutoff
    %if multiple size distributions are input, hix will be a vector...
    %    with length equal to the number of size distributions
%hix2 = highest bin number below upper cutoff
    %for example, if hix2 = 95, there are 95 bins below the upper cutoff
    %inclusive of the bins below the lower cutoff
%mp3 = fit parameters for trimodal distributions
    %number of trimodal size distributions by 9
    %[N0_small, mu_small, lambda_small, N0_mid, mu_mid, lambda_mid,...
    %    N0_large, mu_large, lambda_large]
%mp4 = fit parameters for distributions with small diameter bimodality
    %number of such distributions by 6
    %[N0_small, mu_small, lambda_small, N0_mid, mu_mid, lambda_mid]
%mp5 = fit parameters for distributions with large diameter bimodality
    %number of such distributions by 6
    %[N0_mid, mu_mid, lambda_mid, N0_large, mu_large, lambda_large]
%mp6 = fit parameters for unimodal distributions
    %number of unimodal distributions by 3
    %{N0, mu, lambda]

%note: code uses N(D) = N0*(lambda*D)^mu*exp(-lambda*D)
%to convert to N(D) = N0*D^mu*exp(-lambda*D), multiply N0 by lambda^mu


%Inputs to Input.m:

%time = time of size distribution (s)
%fit_moments = fitting moments used
%bins = bins of size distribution (cm)
%bins_diff = bin widths (cm)
%sd = size distribution (cm^-3 um^-1)
%iwc = ice water content (g m^-3)
%rawcount = raw counts in each bin
%sqError = squared error in bin_concentration (cm^-6)
    %bin_concentration = sd(bins)*bins_diff(bins)
    %note: ensure correct order of magnitude given units above
        %bin_concentration = 10^4*sd(bins)*bins_diff(bins)
%idx = indices (one for each SD) separating smaller habits (sphere, column, 
 %plate) from larger habits (graupel, dendrite, aggregate)
    %if habits are unknown, or if the use of habits to control bimodal2
     %parameters is undesired, input anything here
    %Note: this feature was added while working with data from IMPACTS and
     %might not reflect the typical habit breakdown of modes from HAIC/HIWC

%Additional inputs to TrimodalityTest.m, set in Input.m:

%uonly = flag to force unimodal distributions
uonly = 0;

%testsH = cutoffs to test for large diameter bimodality, in ascending order
testsH = logspace(-1.3, -0.5, 9);
testsH = testsH+1e-5;

%testsL = cutoffs to test for small diameter bimodality, in descending order
testsL = logspace(-1.6, -2.1, 6);
testsL = testsL+1e-5;

%alpha = "sensitivity"
    %suppose you have a trimodal distribution with cutoffs at 0.015 and 0.1 cm
    %alpha is the minimum fraction of N(D)
    %    from the small mode between Dmin and 0.015 cm
    %    from the medium mode between 0.015 and 0.1 cm
    %    from the large mode between 0.1 cm and Dmax
    %sqrt(1-alpha) is the maximum fraction of the expected value
    %    for the ratio between moments 1 and 0 between Dmin and 0.015 cm
    %    for the ratio between moments 2 and 3 between 0.1 cm and Dmax 
alpha = 0.4;

%intmethod = integration method used by TrimodalityTest.m
    %1 = rectangular
    %2 = trapezoidal
    %3 = Simpson's method
    %4 = Simpson's method 3/8
    %5 = Simpson's method composite
intmethod = 1;

%decider = method of choosing cutoffs
    %0 = maximize depth of dip near cutoff
    %1 = minimize overlap between modes
decider = 0;

%conf = confidence level
conf = 0.95;

%cts = minimum # counts between largest small test cutoff and smallest large test cutoff
cts = chi2inv(conf,3);

%lamlam = typical ratio between mid-range and large diameter lambda values
    %used to set bounds on diameter intervals in dip test
    %lamlam = (upper_bound - cutoff)/(cutoff - lower_bound) = (upper_bound - Dmin)/(lower_bound - Dmin)
        %equation in line above used to define upper_bound and lower_bound
lamlam = 2.5; %based on typical lambda values from Darwin

%hdb2 = method controlling whether testsH is arbitrary or derived from habits
    %0 for arbitrary testsH
    %1 for testsH derived from habits (Important: use only if you are confident in habit classification, and be sure to test against hdb2=0)
hdb2 = 0;

if hdb2
    testsH = idx;
end

[tr, ba, bb, un, hix, hix2, mp3, mp4, mp5, mp6] =...
    TrimodalityTest(uonly, testsH, testsL, alpha, time, intmethod, fit_moments,...
    bins, bins_diff, sd, iwc, rawcount, sqError, decider, conf, cts, lamlam, hdb2);
