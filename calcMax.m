function [findDiff] = calcMax(dist, bins, bins_diff)
%remove outlier D in case of bad processing

below = 1;
Dmax = bins(1)+0.5*bins_diff(1);
DmaxX = bins(end)-0.5*bins_diff(end);
while below
    lowPct = bimTest(Dmax, dist, bins, bins_diff, 0);
    if lowPct > erf(2*sqrt(2))
        below = 0;
    end
    if Dmax > DmaxX
        below = 0;
    end
    findDiff = Dmax+0.5*bins_diff-bins;
    [~,findDiff] = min(abs(findDiff));
    %findDiff = find(abs(findDiff) < 0.0001);
    Dmax = Dmax+bins_diff(findDiff);
end
Dmax = min(1.5*Dmax,bins(end)-0.5*bins_diff(end));
findDiff = Dmax+0.5*bins_diff-bins;
[~,findDiff] = min(abs(findDiff));
%findDiff = find(abs(findDiff) <= 0.0051);
%findDiff = max(findDiff);