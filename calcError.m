function error = calcError(sqError, bins, ii)

sz = size(sqError);
error = zeros(sz(1),1);
for j = 1:sz(1)
    error(j) = sqrt(sum(sqError(j,:).*bins.^(2*ii)));
end

end
