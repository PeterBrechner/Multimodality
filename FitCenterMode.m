function [minparams] = FitCenterMode(idx1, idx2, intmethod, sd, bins, bins_diff, lows, error, upr)

%Output: minparams denotes the three fit parameters for the center mode of the distribution

%Fitting to moments 0-2 reduces uncertainty in the large diameter cutoff
%    by reducing uncertainty in the center mode's lambda value
%Use chosen fit moments when testing for small mode
fit_moments = [0 1 2];
if lows
    fit_moments = [0 2 4];
end
sz = size(sd);

minparams = zeros(sz(1),3);

%Fit only bins with diameters less than 0.64 cm
%    Needs to become a non-arbitrary outlier filtering process
%upr = max(find(bins < 0.64));
%sd = sd(:,1:upr);
%bins = bins(1:upr);
%bins_diff = bins_diff(1:upr);
%error = error(:,1:upr);
%idx1 = min(idx1,upr);
%idx2 = min(idx2,upr);

for j=1:sz(1)
    sd(j,upr(j):end) = 0;
    error(j,upr(j):end) = 0;
end

%Calculate observed incomplete moments
M = zeros(sz(1),7);
for ii = 0:6
    for j=1:sz(1)
        M(j,ii+1) = intMethods(intmethod, sd(j,(idx1(j)+1):idx2(j)), bins((idx1(j)+1):idx2(j)), bins_diff((idx1(j)+1):idx2(j)), ii);
    end
end

%Calculate standard deviations for moments 0-2
sigma = zeros(sz(1),3);
for ii = 1:3
    sigma(:,ii) = sigmas(error,bins,idx1,idx2,2,fit_moments(ii));
end

%Fit center modes of size distributions
starting = [2e-3 3 120]; % initial guess
upper = [10 7 280];
lower = [0 -1 0];
%options = optimset('tolfun',1e-16,'tolx',1e-10,'MaxFunEvals',100,'MaxIter',20);
options = optimset('tolfun',1e-4,'tolx',1e-2,'MaxFunEvals',150,'MaxIter',30,'display','off');
for j=1:sz(1)
    if mod(j,1000) == 1
        disp(strcat("FitCenterMode ",num2str(round(100*(j-1)/sz(1))),"% complete: ", string(datetime("now"))));
    end
    %Correct starting first moments for fitting accuracy and precision 
    starting(1) = M(j,2);
    Dcoff1 = bins(1+idx1(j))-0.5*bins_diff(1+idx1(j));
    Dcoff2 = bins(idx2(j))+0.5*bins_diff(idx2(j));
    Dmin = bins(1)-0.5*bins_diff(1);
    Dmax = bins(end)+0.5*bins_diff(end);
    try
        [minparams(j,:), ~, ~, ~, ~] = lsqnonlin(@fit_1_mode, ...
            starting, lower, upper, options, Dmin, Dmax, Dcoff1, Dcoff2, ...
            M(j,fit_moments+1), fit_moments, sigma(j,:)); % force -1 < mu < 5
    catch exception
        disp(sd(j,:));
        rethrow(exception);
    end
    %Calculate N0 from first moment, mu, and lambda
    minparams(j,1) = minparams(j,1)/gamma(minparams(j,2)+2)*minparams(j,3)^(minparams(j,2)+2);
end
