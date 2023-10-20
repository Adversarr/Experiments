clear;
p = readmatrix("minp2w_10_60_300_pr.csv");
r = readmatrix("minp2w_10_60_300_r.csv");
z = readmatrix("minp2w_10_60_300_z.csv");
hz = z(2) - z(1);
hr = r(2) - r(1);
% Gradient in each direction:
%  Because p is 201x100, r is 100x1, z is 201x1. First dimension of P is z. 
[gz, gr] = gradient(p, hz, hr);

% Matlab's `del2` function takes a coefficient 1/(2N). see del2 document.
laplace = 4 * del2(p, hz, hr); 

% Helmholtz equation in R Z coordinate is:
% laplace(P) + 1/r P_r + k^2 P = 0
k = 2 * pi * 300 / 1500;
result = laplace + repmat(1 ./ r', len(z), 1) .* gr;
heatmap(p);
figure;
heatmap(result + k^2 * p);
