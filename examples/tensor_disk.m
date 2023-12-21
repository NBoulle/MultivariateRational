%% Compute multivariate rational approximation on the square
% Reproduce numerical examples in Section 2.4 of the paper

f = @(r,t) cos(10*r+10*pi*t).*(r<=0.75)-sqrt(1-r).*cos(10*r-10*pi*t).*(r>0.75);
g = square_approximation(f, [0,1,-1,1], 'trig_y', 'pole_x', [3/4,1]);

% Evaluate the functions at the grid
Nx = 1000;
Ny = 1000;
r = linspace(0,1,Nx)';
t = [trigpts(Ny-1);1]';
XX = r.*cos(pi*t);
YY = r.*sin(pi*t);
[R,T] = ndgrid(r,t);
F = f(R, T);
F = reshape(F, [Nx, Ny]);
F = F*1.0;
G = g(r,t');
G = reshape(G, [Nx, Ny]);

% Plot f
subplot(1,2,1)
surf(XX,YY,F, 'EdgeColor','interp', 'FaceColor','interp')
colorbar
clim([-1,1])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis square
view(0,90)

% Plot error
subplot(1,2,2)
surf(XX,YY,log10(abs(real(F)-real(G))+eps), 'EdgeColor','interp', 'FaceColor','interp')
colorbar
clim([-16,-12])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis square
view(0,90)
title(sprintf("Error = %.2e", max(abs(F-G),[],"all")))