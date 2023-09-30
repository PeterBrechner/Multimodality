function [ratio, minparams, overlap] = testCutoff(alpha, b, minparams,...
    minparams2, coff, border, inner, lows, intmethod, sd, bins,...
    bins_diff, sqError, iwc, rawcount, conf, cts, lamlam, upr)

%Outputs:
%ratio is the output of the dip test
%minparams is the estimated fit parameters for the center mode
%overlap is the output of the overlap test

%Fit only bins with diameters less than 0.64 cm
%    Needs to become a non-arbitrary outlier filtering process
%upr = max(find(bins < 0.64));
%b = min(b,upr);
%sd = sd(:,1:upr);
%bins = bins(1:upr);
%bins_diff = bins_diff(1:upr);
%rawcount = rawcount(:,1:upr);
%sqError = sqError(:,1:upr);

%Initialize arrays
sz = size(sd);
ratio = ones(sz(1),1);
ratio = (10+coff).*ratio;
Dmin = (bins(1)-0.5*bins_diff(1))*ones(sz(1),1);
Dmax = (bins(end)+0.5*bins_diff(end))*ones(sz(1),1);
for j=1:sz(1)
    sd(j,upr(j):end) = 0;
    rawcount(j,upr(j):end) = 0;
    sqError(j,upr(j):end) = 0;
end
upr
overlap = ones(sz(1),1);
tr = ones(sz(1),1);

%hi is the bin index of the test cutoff
hi = max(find(bins < coff));
vecHi = hi*ones(sz(1),1);

%idx0 is the bin indices between the largest small cutoff...
%    and the smallest large cutoff
idx0 = find(bins > inner(1) & bins < inner(2));

%idxa is the bin indices below the cutoff
idxa = find(bins < coff);

%idxb is the bin indices above the cutoff
idxb = find(bins > coff);

%sums checks for sufficient counts in each mode
sums = zeros(sz(1),1);
for j=1:sz(1)
    sums(j) = sum(rawcount(j,idx0));
    xxxx = find(rawcount(j,idx0));
    if length(xxxx) < 3
        sums(j) = 0;
    end
    if b(j) == sz(2) %no large diameter mode
        lo = min(find(bins > inner(1)));
        sums(j) = sum(rawcount(j,lo:end));
        xxxx = find(rawcount(j,lo:end));
        if length(xxxx) < 3
            sums(j) = 0;
        end
    end
    xxxy = find(rawcount(j,idxa));
    xxxz = find(rawcount(j,idxb));
    if length(xxxy) < 3 | length(xxxz) < 3
        sums(j) = 0;
    end
end

%Redefining idx0 to be the distributions with sufficient counts
idx0 = find(sums >= cts);

%Estimate gamma fit parameters for the center mode
if coff < border
%lambda for the small diameter mode is inversely proportional to coff
%    so we can take advantage
lamlam = lamlam*inner(1)/coff;
minparams(idx0,:) = FitCenterMode(alpha, vecHi(idx0), b(idx0), intmethod,...
    sd(idx0,:), bins, bins_diff, lows, sqError(idx0,:), upr);
end
if coff > border
minparams(idx0,:) = FitCenterMode(alpha, b(idx0), vecHi(idx0), intmethod,...
    sd(idx0,:), bins, bins_diff, lows, sqError(idx0,:), upr);
end

%Store gamma fit parameters for center mode as N0, mu, lambda
N0 = minparams(:,1);
mu = minparams(:,2);
lambda = minparams(:,3);

%Store gamma fit parameters for dip test as N01, mu1, lambda1
N01 = minparams2(:,1);
mu1 = minparams2(:,2);
lambda1 = minparams2(:,3);

%Perform tests
for j=1:sz(1)
    %Only when counts are sufficient
    if sums(j) < cts
        continue
    end
    %Overlap test
    if coff < border
        Dmax(j) = bins(b(j))+0.5*bins_diff(b(j));
        overlap(j) = overlapTest(Dmin(j), coff, N0(j), mu(j), lambda(j),...
            alpha, sd(j,:), bins, bins_diff,...
            sqError(j,:), conf, 0, 0);
    end
    if coff > border
        Dmin(j) = bins(b(j)+1)-0.5*bins_diff(b(j)+1);
        overlap(j) = overlapTest(coff, Dmax(j), N0(j), mu(j), lambda(j),...
            alpha, sd(j,:), bins, bins_diff,...
            sqError(j,:), conf, 0, 2);
    end
    if overlap(j) >= alpha
        tr(j) = 0;
    end
    %Dip test
    coff2 = 0.5*coff-0.5*Dmin(j);
    idx1 = find(bins > Dmin(j) & bins < coff-coff2*(1-1/lamlam));
    idx2 = find(bins >= coff-coff2*(1-1/lamlam) & bins <= coff+coff2*(lamlam-1));
    idx3 = find(bins > coff+coff2*(lamlam-1) & bins < Dmax(j));
    ll = bins(idx1(1))-0.5*bins_diff(idx1(1));
    ml = bins(idx2(1))-0.5*bins_diff(idx2(1));
    ul = bins(idx3(1))-0.5*bins_diff(idx3(1));
    uu = bins(idx3(end))+0.5*bins_diff(idx3(end));
    if tr(j) & sum(rawcount(j,idx1)) >= 0.5 & sum(rawcount(j,idx2)) >= 0 & sum(rawcount(j,idx3)) >= 0.5
        [ratio(j),overlap(j)] = dipTest(idx1, idx2, idx3, ll, ml, ul, uu,...
            N01(j), mu1(j), lambda1(j), conf, alpha, overlap(j), sd(j,:),...
            bins, bins_diff, sqError(j,:), 0);
    end
end

end
