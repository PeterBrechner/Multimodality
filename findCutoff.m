function [coff, minr, b, mp] = findCutoff(alpha, tests, border, b, mp, mp2,...
    lows, intmethod, sd, bins, bins_diff, sqError, rawcount, decider,...
    conf, cts, lamlam, inner, upr)

%Outputs:
%coff = cutoff diameter
%minr = deepest dip or minimum overlap
%b = bin index for cutoff diameter
%mp = fit parameters for center mode

%Inputs:
%b = bin index for opposite cutoff diameter

%Initialize arrays
sz = size(sd);
coff = zeros(sz(1),1);
minr = coff;
ratios = zeros(sz(1),length(tests));
overlap = zeros(sz(1),length(tests));
mps = zeros(sz(1),3*length(tests));
inner = [min(inner, tests(1)), max(inner, tests(1))];

%Test cutoffs
for c = 1:length(tests)
    disp(strcat("Testing Dcutoff ",num2str(tests(c))," cm"));
    [ratios(:,c), mps(:,(3*c-2):3*c), overlap(:,c)] =...
        testCutoff(alpha, b, mp, mp2, tests(c), border, inner, lows,...
        intmethod, sd, bins, bins_diff, sqError, rawcount, conf, cts,...
        lamlam, upr);
    if decider
        ratios(:,c) = -alpha+overlap(:,c);
    end
end

%ratios
%overlap

%Determine cutoffs
for j=1:sz(1)
    coff(j) = tests(find(ratios(j,:) == min(ratios(j,:)), 1));
    minr(j) = min(ratios(j,:));
    mp(j,:) = mps(j,(3*find(ratios(j,:) == minr(j), 1)-2):3*find(ratios(j,:) == minr(j), 1));
    if minr(j) < 0
        b(j) = find(bins < coff(j),1,'last');
    else
        if coff(j) < border
            b(j) = 0;
        else
            b(j) = sz(2);
        end
        mp(j,:) = mp2(j,:);
    end
end
