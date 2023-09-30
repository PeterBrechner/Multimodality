function lowPct = bimTest(coff, sd, bins, bins_diff, ii)

idx1 = find(bins < coff);
idx2 = find(bins > coff);

lower = intMethods(1, sd(:,idx1), bins(idx1), bins_diff(idx1), ii);
upper = intMethods(1, sd(:,idx2), bins(idx2), bins_diff(idx2), ii);

lowPct = lower./(lower+upper);

end

