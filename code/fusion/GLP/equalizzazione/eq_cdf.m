function v1 = eq_cdf( v1, v2, N_pts )
[val1, cdf1] = calc_cdf (v1, N_pts);
[val2, cdf2] = calc_cdf (v2, N_pts);

for i=1:length(cdf1)
    dist = abs( cdf2 - cdf1(i));
    [tmp, j] = min (dist);
    if i==1
        v1(v1<=val1(i)) = val2(j);
    else
        v1(v1>val1(i-1) & v1<=val1(i)) = val2(j);
    end
end
end

