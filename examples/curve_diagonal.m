%% Example of a diagonal singularity
%
% The following code is an example of the approximation of a function with
% a singularity along the diagonal. Here, we use a polynomial basis in x+y
% for the varying residue.
% This example is included for information purposes and is not a part of
% the paper. The methodology is described in section 4.

f = @(x,y) sqrt(abs(x-y)) .* exp(x+y);

%% Cluster points in the normal direction of the characteristic curves.

threshold = 1e-12;

% Define poles clustering towards zero
Nq = 100;
sigma = 4;
qj = exp(-sigma*(1:Nq)/sqrt(Nq));

% Define domain
dom_x = [0,1];
dom_y = [0,1];

% Curve
f_curve = @(x,y) x-y;

Ns = 15;        % Degree of smooth part in the x and y coordinates
Np = 15;         % Degree of smooth part for the varying residue

% Construct 2d Chebyshev grid
Ms = 3*Ns;
X1 = chebpts(Ms, dom_x);
Y1 = chebpts(Ms, dom_y);
[x1,y1] = ndgrid(X1,Y1);

% Constructs points clustering towards the diagonal
Mq = 3*Nq;       % nb of points in the normal direction
Mp = 3*Np;       % nb of points in the tangent direction
[x2,y2] = compute_clustered_points(f_curve, Mp, Mq, [dom_x,dom_y], threshold);

% put all of the points together in a long vector
X = [x1(:); x2(:)];
Y = [y1(:); y2(:)];


%% Evaluate the function at the grid of points

% We also make a plotting grid
tx = linspace(dom_x(1), dom_x(2), 1000)';
ty = linspace(dom_y(1), dom_y(2), 1000)';
[Tx,Ty] = meshgrid(tx, ty);

F = f(X, Y);
Fplot = f(Tx,Ty);


%% Now continue with the approximation

% Get + and - imaginary poles
qj_x = [];
for q = qj
    qj_x = [qj_x, 1i*q, -1i*q];
end

% Define basis of rational functions
r = @(x,y) qj_x ./ (x - y + eps + qj_x);             % partial fractions
Psi = r(X,Y);
pxy = chebpoly(0:Np-1,[0, dom_x(end) + dom_y(end)]); % varying residue
Phi_xy = pxy(X + Y);
Psi_poly = linearize_tensorproduct(Psi, Phi_xy);

% Define 2d tensor product Chebyshev polynomials
px = chebpoly(0:Ns-1,dom_x);
py = chebpoly(0:Ns-1,dom_y);
Phi_x = px(X);
Phi_y = py(Y);
Phi_poly = linearize_tensorproduct(Phi_x, Phi_y);

% Setup and solve the system
A = [Phi_poly Psi_poly];
c = A\F;
fprintf("Residual: %d\n", norm(A*c-F))
fprintf("Coefficient norm: %d\n", norm(c))

%% Make functions to evaluate the solution at a point
c_poly = reshape(c(1:Ns^2),Ns,Ns);
c_rat = reshape(c(Ns^2+1:end),2*Nq,Np);

% these evaluation routines take matrices X and Y of equal size
cp = @(X,Y) arrayfun(@(x,y) px(x) * c_poly * py(y)', X, Y);
cr = @(X,Y) arrayfun(@(x,y) r(x,y) * c_rat * pxy(x+y)', X, Y);

cf = @(x,y) cp(x,y)+cr(x,y);

% the polynomials have product structure so we can also work with vectors
cp_vec = @(xvec,yvec) px(xvec) * c_poly * py(yvec)';

Zpoly = cp_vec(tx,ty).';
Zrat = cr(Tx,Ty);
Zplot = Zpoly + Zrat;

%% Plot the samples

figure
subplot(2,2,1)
plot(X, Y, 'o')
axis equal, xlabel("x"), ylabel("y")
title("sample points")
xlim(dom_x); ylim(dom_y);

subplot(2,2,2)
plot(x1(:), y1(:), 'o')
axis equal, xlabel("x"), ylabel("y")
title("Chebyshev points")
xlim(dom_x); ylim(dom_y);

subplot(2,2,3)
plot(x2(:), y2(:), 'o')
axis equal, xlabel("x"), ylabel("y")
title("clustered points")
xlim(dom_x); ylim(dom_y);

subplot(2,2,4)
plot3(X, Y, F, 'o')
axis equal, xlabel("x"), ylabel("y")
title("function samples")

%% Plot the approximation


figure
subplot(3,2,1)
surf(tx, ty, Fplot)
shading interp
title("function to approximate")

% the approximation
subplot(3,2,2)
surf(tx, tx, real(Zplot))
shading interp
title("the approximation")

subplot(3,2,3)
surf(tx, ty, real(Zpoly))
shading interp
title("the polynomial part")

subplot(3,2,4)
surf(tx, ty, real(Zrat))
shading interp
title("the rational part")

% the approximation error
subplot(3,2,5)
surf(tx, ty, log10(abs(Fplot-Zplot)))
title("approximation error")
shading interp

shg

