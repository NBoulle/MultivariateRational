%% Piecewise multivariate rational approximation
% Reproduce numerical example in Section 3
% Require the ultraSEM package at https://github.com/danfortunato/ultraSEM

% Define function to be approximated
u = @(x,y) cos(5*pi*(x+y)).*sqrt(abs(x-y));

% Define domain
v1 = [0 0 ; 1 1 ; 1 0];
dom1 = ultraSEM.triangle(v1);
v2 = [0 0 ; 1 1 ; 0 1];
dom2 = ultraSEM.triangle(v2);
dom = merge([dom1, dom2]);

% Extract number of patches
d = dom.domain;
kk = length(dom);

% Records points on physical domain
x = cell(kk);
y = cell(kk);
u_eval = cell(kk);
u_true = cell(kk);

% Create points to plot the function on mapped domain
nplotpts = 500;
x1 = linspace(-1, 1, nplotpts);
[xk, yk] = meshgrid(x1);

% Loop over the number of patches
for k = 1:kk
    
    sprintf("k = %d / %d", k, kk)

    % Extract patch k
    d_k = d(k,:);

    % Transform points to physical domain
    [x{k,1}, y{k,1}] = transformGrid(d_k, xk, yk);

    % Map the unit domain to function on the patch
    f = @(r,s) u(d_k.x(r,s), d_k.y(r,s));
    
    % Compute approximation on [-1,1,-1,1]
    if k == 1 || k==3
        g = square_approximation(f, [-1,1,-1,1], 'pole_x', [-1], 'N2', 25, 'N3', 25);
    elseif k==4 || k==5
        g = square_approximation(f, [-1,1,-1,1], 'pole_y', [-1], 'N2', 25, 'N3', 25);
    else
        g = square_approximation(f, [-1,1,-1,1], 'N2', 25, 'N3', 25);
    end

    % Evaluate the approximation on the unit points
    u_eval{k} = g(x1',x1')'; 
    u_true{k} = f(xk, yk);
end

% Plot the approximation on the physical domain
close all
subplot(1,2,1)

for k = 1:kk
    surf(x{k}, y{k}, real(u_eval{k}), 'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on
end
hold off
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
colorbar
clim([-1,1])
axis square
view(0,90)

% Plot errors

subplot(1,2,2)
for k = 1:kk
    surf(x{k}, y{k}, log10(abs(u_true{k}-real(u_eval{k}))+eps), 'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on
    max(max(abs(u_true{k}-real(u_eval{k}))))
end
hold off
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
colorbar
axis square
view(0,90)