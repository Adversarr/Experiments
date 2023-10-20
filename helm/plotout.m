clear;
p = readmatrix("output.txt");
% heatmap(m); 
R = 50;
Z = 50;
[nr, nz] = size(p);
hr = R / (nr - 1);
hz = Z / (nz - 1);
r = linspace(0, R, nr);
z = linspace(0, Z, nz);
[gz, gr] = gradient(p, hz, hr);
% Matlab's `del2` function takes a coefficient 1/(2N). see del2 document.
laplace = 4 * del2(p, hz, hr); 
% Helmholtz equation in R Z coordinate is:
% laplace(P) + 1/r P_r + k^2 P = 0
k = 2 * pi * 50 / 1500;
result = laplace + repmat(1 ./ (r + 1e-7), nz, 1) .* gr;
% [R, Z] = meshgrid(r, z);
imagesc(r, z, p);
title("p");
colorbar;
figure;
err = (result + k^2 * p);
heatmap(r(2:end-1), z(2:end-1), log(abs(err(2:end-1, 2:end-1))));
colorbar;
title("Log |\Delta p + k^2 p|");


