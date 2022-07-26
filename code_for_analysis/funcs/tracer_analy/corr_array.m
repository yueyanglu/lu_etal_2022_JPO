function r = corr_array(a,b,dim)

%--- Compute correlations on third dimension
% Remove means 
az = bsxfun(@minus, a, mean(a,dim));
bz = bsxfun(@minus, b, mean(b,dim));
% variance
a2 = az .^ 2;
b2 = bz .^ 2;
ab = az .* bz;
% Standard Pearson correlation coefficient formula
r = sum(ab, dim) ./ sqrt(sum(a2, dim) .* sum(b2, dim));

end