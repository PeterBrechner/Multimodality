function [minparams] = FitCenterMode(alpha, idx1, idx2, intmethod, sd, bins, bins_diff, fit_moments, error, upr)

sz = size(sd);

minparams = zeros(sz(1),3);
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

M = zeros(sz(1),7);
for ii = 0:6
    for j=1:sz(1)
        M(j,ii+1) = intMethods(intmethod, sd(j,(idx1(j)+1):idx2(j)), bins((idx1(j)+1):idx2(j)), bins_diff((idx1(j)+1):idx2(j)), ii);
    end
end

sigma = zeros(sz(1),3);
for ii = 1:3
    sigma(:,ii) = sigmas(error,bins,idx1,idx2,2,fit_moments(ii));
end

starting = [2e-3 3 120]; % initial guess
upper = [10 9 400];
%lower = [1e-6 -1 0];
lower = [0 -1 0];
options = optimset('tolfun',1e-16,'tolx',1e-10,'MaxFunEvals',100,'MaxIter',20);
for j=1:sz(1)
    j
    starting(1) = M(j,2);
    Dcoff1 = bins(1+idx1(j))-0.5*bins_diff(1+idx1(j));
    Dcoff2 = bins(idx2(j))+0.5*bins_diff(idx2(j));
    Dmin = bins(1)-0.5*bins_diff(1);
    Dmax = bins(end)+0.5*bins_diff(end);
    [minparams(j,:), minchisq, ~, exitflag, ~] =...
        lsqnonlin(@fit_1_mode, starting, lower,...
        upper, options, Dmin, Dmax, Dcoff1, Dcoff2, M(j,fit_moments+1),...
        fit_moments, sigma(j,:)); % force -1 < mu < 5
    minparams(j,1) = minparams(j,1)/gamma(minparams(j,2)+2)*minparams(j,3)^(minparams(j,2)+2);
    ginc0 = gammainc(minparams(j,3)*Dmax,minparams(j,2)+1)-gammainc(minparams(j,3)*Dmin,minparams(j,2)+1);
    g0 = ginc0*minparams(j,1)*gamma(minparams(j,2)+1)/minparams(j,3)^(minparams(j,2)+1)
    a0 = M(j,1)
    ginc2 = gammainc(minparams(j,3)*Dmax,minparams(j,2)+3)-gammainc(minparams(j,3)*Dmin,minparams(j,2)+3);
    g2 = ginc2*minparams(j,1)*gamma(minparams(j,2)+3)/minparams(j,3)^(minparams(j,2)+3)
    a2 = M(j,3)
    ginc4 = gammainc(minparams(j,3)*Dmax,minparams(j,2)+5)-gammainc(minparams(j,3)*Dmin,minparams(j,2)+5);
    g4 = ginc4*minparams(j,1)*gamma(minparams(j,2)+5)/minparams(j,3)^(minparams(j,2)+5)
    a4 = M(j,5)
end
