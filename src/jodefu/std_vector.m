function x = std_vector(x, dim)

dim_total = 1 : ndims(x);
dim_excluded = dim_total;
dim_excluded(dim) = [];
dim_concatenate = cat(2, dim, dim_excluded);
x_permute = permute(x, dim_concatenate);
sizes = size(x);
dim_reshape = cat(2, prod(sizes(dim)), sizes(dim_excluded));
x_reshaped = reshape(x_permute, dim_reshape);
x_std = std(x_reshaped, 0, 1);
sizes(dim) = 1;
x = reshape(shiftdim(x_std, -1), sizes);