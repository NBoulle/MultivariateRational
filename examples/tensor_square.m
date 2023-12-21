%% Compute multivariate rational approximation on the square
% Reproduce numerical examples in Section 2.4 of the paper

% Example 1:
f = @(x,y) (x .* (1-x)).^(1/4+y).*sqrt(y.*(1-y));
g = square_approximation(f, [0,1,0,1], 'pole_x', [0,1], 'pole_y', [0,1]);

% Plot results
X = linspace(0,1,1000)';
Y = linspace(0,1,1000)';
[x,y] = ndgrid(X,Y);
F = f(X,Y');
G = g(X,Y);

close all
subplot(1,2,1)
surf(x,y,real(F), 'EdgeColor','interp', 'FaceColor','interp')
colorbar
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
clim([0,0.25])
axis square
view(0,90)

subplot(1,2,2)
surf(x,y,log10(abs(real(F)-real(G))+eps), 'EdgeColor','interp', 'FaceColor','interp')
colorbar
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
c = colorbar;
clim([-16,-12])
axis square
set(gca,'TickLabelInterpreter','latex')
title(sprintf("Error = %.2e", max(abs(F-G),[],"all")))
view(0,90)

%%
% Example 2
f = @(x,y)sqrt(x+y);
g = square_approximation(f, [0,1,0,1], 'pole_x', [0], 'pole_y', [0]);

% Plot results
X = linspace(0,1,1000)';
Y = linspace(0,1,1000)';
[x,y] = ndgrid(X,Y);
F = f(X,Y');
G = g(X,Y);

close all
subplot(1,2,1)
surf(x,y,real(F), 'EdgeColor','interp', 'FaceColor','interp')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
colorbar
clim([0,1.5])
axis square
view(0,90)

subplot(1,2,2)
surf(x,y,log10(abs(real(F)-real(G))+eps), 'EdgeColor','interp', 'FaceColor','interp')
colorbar
clim([-16,-12])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis square
title(sprintf("Error = %.2e", max(abs(F-G),[],"all")))
view(0,90)