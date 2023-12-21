%% Convergence of the rational approximant
% Reproduce Fig. 3
f = @(x,y)sqrt(x+y);

% Plot results
X = linspace(0,1,1000)';
Y = linspace(0,1,1000)';
[x,y] = ndgrid(X,Y);
F = f(X,Y');

%% Panel (a)
subplot(1,2,1)
ListN1 = 0:5:225;
N2 = 16;
Error = [];
for N1=ListN1
    [g, poles, C] = square_approximation(f, [0,1,0,1], 'pole_x', [0], 'pole_y', [0], "N1", N1, "N2", N2, "N3", N2, "M", 300);
    G = g(X,Y);
    E = max(abs(F-G),[],"all");
    Error = [Error, E];
end
xlabel('$\sqrt{N_1}$','interpreter','latex')
ylabel('$\|f-r_{N_1}\|_{\max}$','interpreter','latex')
semilogy(sqrt(ListN1),Error,'.-')

axis square
imin=2;
imax=14;
a = polyfit(sqrt(ListN1(imin:imax)), log10(Error(imin:imax)),1);
hold on
semilogy(sqrt(ListN1), 10.^(a(1)*sqrt(ListN1)+a(2)))
hold off
ylim([1e-15,1])

%% Panel (b)
ListN1 = 0:5:150;
ListN2 = 0:2:12;
A = [];
for N2 = ListN2
    Error = [];
    for N1=ListN1
        [g, poles, C] = square_approximation(f, [0,1,0,1], 'pole_x', [0], 'pole_y', [0], "N1", N1, "N2", N2, "N3", N2, "M", 300);
        G = g(X,Y);
        E = max(abs(F-G),[],"all");
        Error = [Error, E];
    end
    A = [A; Error];
end

% Contour plot
subplot(1,2,2)
contourf(sqrt(ListN1),ListN2,log10(A)); colorbar
hold on
plot(sqrt(ListN1),1.3*sqrt(ListN1),'w--',"LineWidth",2)
hold off
xlim([0,12])
ylim([0,12])
clim([-15,0])
xlabel('$\sqrt{N_1}$','interpreter','latex')
ylabel('$N_2$','interpreter','latex')
axis square
set(gca,'TickLabelInterpreter','latex')
view(0,90)