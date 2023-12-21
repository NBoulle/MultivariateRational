function A = evaluate_approximant(X, Y, N1, N2, N3, dom_x, dom_y, f_curve)
% Evaluate the rational approximant at points X,Y

% Define poles
sigma = 3;
pj = exp(-sigma*(sqrt(N1) - sqrt(1:N1)));
pj_imag = [1i*pj, -1i*pj];
pj = [pj, pj];

% Define basis of rational functions
r = @(x,y) pj ./ (f_curve(x,y) + eps + pj_imag);
Psi = r(X,Y);

% Define 2d tensor product Chebyshev polynomials
px = chebpoly(0:N2-1,dom_x);
py = chebpoly(0:N2-1,dom_y);
Phi_x = px(X);
Phi_y = py(Y);

Phi_poly = zeros(length(X), length(px)*length(py));
for i = 1:length(py)
    % we linearize assuming x varies along rows and y along columns
    Phi_poly(:,(i-1)*N2+1:i*N2) = Phi_x .* repmat(Phi_y(:,i), 1, N2);
end

% Define 2d tensor product Chebyshev polynomials
px = chebpoly(0:N3-1,dom_x);
py = chebpoly(0:N3-1,dom_y);
Phi_x = px(X);
Phi_y = py(Y);

Phi_poly2 = zeros(length(X), length(px)*length(py));
for i = 1:length(py)
    % we linearize assuming x varies along rows and y along columns
    Phi_poly2(:,(i-1)*N3+1:i*N3) = Phi_x .* repmat(Phi_y(:,i), 1, N3);
 end

Psi_poly = [];
for i = 1:size(Phi_poly2,2)
    Psi_poly = [Psi_poly, Psi .* repmat(Phi_poly2(:,i), 1, length(pj))];
end

A = [Phi_poly Psi_poly];
end