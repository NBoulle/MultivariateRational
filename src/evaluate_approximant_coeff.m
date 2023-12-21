function [Phi_poly,Psi_poly]  = evaluate_approximant_coeff(X, Y, N1, N2, N3, dom_x, dom_y, f_curve, c)
% Given the coefficient vector c, return the corresponding polynomial and
% rational parts evaluated at the points (X,Y).

% Define poles
sigma = 3;
pj = exp(-sigma*(sqrt(N1) - sqrt(1:N1)));
pj_imag = [1i*pj, -1i*pj];
pj = [pj, pj];

% Define basis of rational functions
R = f_curve(X,Y);

% Define 2d tensor product Chebyshev polynomials
px = chebpoly(0:N2-1,dom_x);
py = chebpoly(0:N2-1,dom_y);
Px = px(X,:);
Py = py(Y,:);

Phi_poly = zeros(length(X),1);
k = 1;
for j = 1:length(py)
    for i = 1:length(px)
        % we linearize assuming x varies along rows and y along columns
        Phi_poly = Phi_poly + Px(:,i).*Py(:,j)*c(k);
        k = k+1;
    end
end


% Define 2d tensor product Chebyshev polynomials
px = chebpoly(0:N3-1,dom_x);
py = chebpoly(0:N3-1,dom_y);
Px = px(X,:);
Py = py(Y,:);

Psi_poly = zeros(length(X),1);
for j = 1:length(py)
    for i = 1:length(px)
        for l = 1:length(pj)
        % we linearize assuming x varies along rows and y along columns
        Psi_poly = Psi_poly + Px(:,i).*Py(:,j).*(pj(l) ./ (R + eps + pj_imag(l)))*c(k);
        k = k+1;
        end
    end
end
end