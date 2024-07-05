%% Multivariate rational approximation with elliptic curve of singularity
% Reproduce numerical example in section 4.2

% Define the domain
dom_x = [-2,2];
dom_y = [-2,2];

% Function f to be approximated
f = @(x,y) abs(x.^3-2*x+1-y.^2);

% Algebraic equation of the singularity curve
f_curve = @(x,y) x.^3-2*x+1-y.^2;

% Size of the approximation space
N1 = 50;        % Degree of rational expansion in normal direction
N2 = 60;        % Degree of smooth part in the x and y coordinates
N3 = 3;         % Degree of smooth part for e varying residue

% Construct 2d Chebyshev grid
n_cheb = N2;
X1 = chebpts(n_cheb, dom_x);
Y1 = chebpts(n_cheb, dom_y);
[x1,y1] = ndgrid(X1,Y1);

n_rho = 2*N1;
n_theta = 20;    % nb of points in the angular direction

% Construct clustering grid of points
[x2, y2] = compute_clustered_points(f_curve, n_theta, n_rho, [dom_x,dom_y]);

% Put all of the points together in a long vector
X = [x1(:); x2(:)];
Y = [y1(:); y2(:)];

% Evaluate the function at the grid of points
F = f(X,Y);

% Evaluate the rational approximant
A = evaluate_approximant(X, Y, N1, N2, N3, dom_x, dom_y, f_curve);

% Setup and solve the system
tic; c = A\F; toc;
disp(sprintf("Residual: %d", norm(A*c-F)))
disp(sprintf("Coefficient norm: %d", norm(c)))

%% Extract rational and smooth residue parts at grid of points
% Sample points
nplot = 1000;
x = linspace(dom_x(1), dom_x(2),nplot);
y = linspace(dom_y(1), dom_y(2),nplot);
[x1plot,y1plot] = ndgrid(x,y);
Xplot = x1plot(:);
Yplot = y1plot(:);
[Phi_poly,Psi_poly] = evaluate_approximant_coeff(Xplot, Yplot, N1, N2, N3, dom_x, dom_y, f_curve, c);

%% Plot the function we're approximating

% Evaluate the approximation
close all
subplot(2,2,1)
surf(x1plot, y1plot, reshape(real(Phi_poly+Psi_poly),nplot,nplot))
shading interp
colorbar
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
clim([0,8])
zlim([0,8])
axis square
set(colorbar,'visible','off')
title("Rational Approximation")

% Plot the error
subplot(2,2,2)
E = abs(reshape(real(Phi_poly+Psi_poly), nplot, nplot)-f(x1plot,y1plot));
surf(x1plot, y1plot, log10(E))
shading interp
colorbar
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
clim([-12,-8])
zlim([-12,-8])
title(sprintf("Error: %.2e", max(E(:))))

% Plot polynomial part
subplot(2,2,3)
surf(x1plot, y1plot, real(reshape(Phi_poly, nplot, nplot)))
shading interp
axis square
view(0,90)
colorbar
clim([1,7])
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
title("Polynomial part")

% Plot rational part
subplot(2,2,4)
surf(x1plot, y1plot, real(reshape(Psi_poly, nplot, nplot)))
shading interp
axis square
view(0,90)
colorbar
clim([-1,0])
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
title("Rational part")
