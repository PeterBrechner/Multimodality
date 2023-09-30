function error = sigmas(squaredErrorCount, bins, coff1, coff2, mode, ii)

sz = size(squaredErrorCount);
error = zeros(sz(1),1);
for j = 1:sz(1)
    if mode == 1
        sEC = squaredErrorCount(j,1:coff1(j));
        b = bins(1:coff1(j));
        error(j) = sqrt(sum(sEC.*b.^(2*ii)));
    elseif mode == 2
        sEC = squaredErrorCount(j,(1+coff1(j)):coff2(j));
        b = bins((1+coff1(j)):coff2(j));
        error(j) = sqrt(sum(sEC.*b.^(2*ii)));
    else
        sEC = squaredErrorCount(j,(1+coff2(j)):end);
        b = bins((1+coff2(j)):end);
        error(j) = sqrt(sum(sEC.*b.^(2*ii)));
    end
end

end
