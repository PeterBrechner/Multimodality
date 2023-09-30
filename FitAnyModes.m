function [minparams] = FitAnyModes(alpha, mag, idx1, idx2, intmethod, sd, bins, bins_diff, fit_moments, error, upr)

sz = size(sd);

minparams = zeros(sz(1),6);
%upr = max(find(bins < 0.64));
%sd = sd(:,1:upr);
%bins = bins(1:upr);
%bins_diff = bins_diff(1:upr);
%error = error(:,1:upr);
%idx1 = min(idx1,upr);
%idx2 = min(idx2,upr);

for j=1:sz(1)
    if idx2(j) == bins(end)+0.5*bins_diff(end)
        upr(j) = sz(2);
    else
        sd(j,upr(j):end) = 0;
        error(j,upr(j):end) = 0;
    end
end

M = zeros(sz(1),7);
M2 = M;
M3 = M;
for ii = 0:6
    for j=1:sz(1)
        if idx1(j) ~= 0
            M(j,ii+1) = intMethods(intmethod, sd(j,1:idx1(j)),...
                bins(1:idx1(j)), bins_diff(1:idx1(j)), ii);
        end
        M2(j,ii+1) = intMethods(intmethod, sd(j,(idx1(j)+1):idx2(j)),...
            bins((idx1(j)+1):idx2(j)), bins_diff((idx1(j)+1):idx2(j)), ii);
        if idx2(j) ~= upr(j)
            M3(j,ii+1) = intMethods(intmethod, sd(j,(idx2(j)+1):end),...
                bins((idx2(j)+1):end), bins_diff((idx2(j)+1):end), ii);
        end
    end
end

sigma = zeros(sz(1),3);
sigma2 = zeros(sz(1),3);
sigma3 = zeros(sz(1),3);
for ii = 1:3
    sigma(:,ii) = sigmas(error,bins,idx1,idx2,1,fit_moments(ii));
    sigma2(:,ii) = sigmas(error,bins,idx1,idx2,2,fit_moments(ii));
    sigma3(:,ii) = sigmas(error,bins,idx1,idx2,3,fit_moments(ii));
end

%options = optimset('tolfun',1e-16,'tolx',1e-10,'MaxFunEvals',1575,'MaxIter',75);
options = optimset('tolfun',1e-4,'tolx',1e-2,'MaxFunEvals',1140,'MaxIter',60,'display','off');
for j=1:sz(1)
    upper = [1, 9, 2000, 10, 9, 400, 1, 9, 400];
    lower = [0, -1, 200, 0, -1, 0, 0, -1, 0];
    Dcoff1 = bins(1+idx1(j))-0.5*bins_diff(1+idx1(j));
    Dcoff2 = bins(idx2(j))+0.5*bins_diff(idx2(j));
    Dmin = bins(1)-0.5*bins_diff(1);
    Dmax = bins(end)+0.5*bins_diff(end);
    starting = find_starting(mag(j,:), Dcoff1, Dcoff2, Dmin, Dmax);
    if idx1(j) == 0
        if mod(j,100) == 1
            disp(strcat("Bimodal2 Fits ",num2str(round(100*(j-1)/sz(1))),"% complete: ", string(datetime("now"))));
        end
        upper = upper(4:9);
        lower = lower(4:9);
        [minparams(j,:), ~, ~, ~, ~] =...
            lsqnonlin(@fit_bimodal2, starting, lower,...
            upper, options, alpha, Dmin, Dmax, Dcoff2,...
            M2(j,1:5), M3(j,1:5),...
            fit_moments, sigma2(j,:), sigma3(j,:));
%        minparams = starting;
        minparams(j,1) = minparams(j,1)/gamma(minparams(j,2)+2)*minparams(j,3)^(minparams(j,2)+2);
        minparams(j,4) = minparams(j,4)/gamma(minparams(j,5)+2)*minparams(j,6)^(minparams(j,5)+2);
    else
        if mod(j,100) == 1
            disp(strcat("Bimodal1 Fits ",num2str(round(100*(j-1)/sz(1))),"% complete: ", string(datetime("now"))));
        end
        upper = upper(1:6);
        lower = lower(1:6);
        [minparams(j,:), ~, ~, ~, ~] =...
            lsqnonlin(@fit_2_modes, starting, lower,...
            upper, options, alpha, Dmin, Dmax, Dcoff1,...
            M(j,1:5), M2(j,1:5),...
            fit_moments, sigma(j,:), sigma2(j,:));
%        minparams = starting;
        minparams(j,1) = minparams(j,1)/gamma(minparams(j,2)+2)*minparams(j,3)^(minparams(j,2)+2);
        minparams(j,4) = minparams(j,4)/gamma(minparams(j,5)+2)*minparams(j,6)^(minparams(j,5)+2);
    end
end
