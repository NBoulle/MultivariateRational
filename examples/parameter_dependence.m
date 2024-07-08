%% Show the dependence of the accuracy on the default parameter settings
% Reproduce numerical examples in Section 2.5

f1 = @(x,y) (x .* (1-x)).^(1/4+y).*sqrt(y.*(1-y));
f2 = @(x,y) sqrt(x+y);
f3 = @(r,t) cos(10*r+10*pi*t).*((r<=0.75)-sqrt(1-r).*(r>0.75));

%% Dependence of accuracy on sigma
% independent grid for error computation
M = 1000; % nb of grid pts per pole ( compared to 5N_1 sample pts)
%   --> f1 (poles at x=0,1 and y=0,1)
X1 = [logspace(log10(eps),0,M)'; 1 - logspace(log10(eps),0,M)'; chebpts(M,[0,1])];
Y1 = [logspace(log10(eps),0,M)'; 1 - logspace(log10(eps),0,M)'; chebpts(M,[0,1])];
F1 = f1(X1,Y1');
%   --> f2 (poles at x=0 and y=0)
X2 = [logspace(log10(eps),0,M)'; chebpts(1000,[0,1])];
Y2 = [logspace(log10(eps),0,M)'; chebpts(1000,[0,1])];
F2 = f2(X2,Y2');
%   --> f3 (poles at r=0.75,1)
r = [0.75 + 0.25*logspace(log10(eps),0,M/2)'; 0.75 - ...
    0.75*logspace(log10(eps),0,M/2)'; 1 - logspace(log10(eps),0,M)'; chebpts(M,[0,1])];
r = r(r ~= 0.75);           % delete chebpt at 0.75
t = [trigpts(M-1); 1]';
[R3,T3] = ndgrid(r,t);
F3 = 1.0*reshape(f3(R3, T3), [3*M-1, M]);

ListSigma = linspace(0,20,300);
Error_1 = []; Error_2 = []; Error_3 = [];
for sigma = ListSigma
    [g, ~, ~] = square_approximation(f1, [0,1,0,1], 'pole_x', [0,1], ...
        'pole_y', [0,1], 'sigma', sigma, 'N1', 40);
    G = g(X1,Y1);
    E = max(abs(F1-G),[],"all");
    Error_1 = [Error_1, E];

    [g, ~, ~] = square_approximation(f2, [0,1,0,1], 'pole_x', [0], ...
        'pole_y', [0], 'sigma', sigma, 'N1', 40);
    G = g(X2,Y2);
    E = max(abs(F2-G),[],"all");
    Error_2 = [Error_2, E];

    [g, ~, ~] = square_approximation(f3, [0,1,-1,1], 'trig_y', ...
        'pole_x', [3/4,1], 'sigma', sigma, 'N1', 80);
    G = g(r,t');
    E = max(abs(F3-G),[],"all");
    Error_3 = [Error_3, E];

    sprintf("sigma = %.2e,", sigma)
end

figure;
semilogy(ListSigma,Error_1,'.-r'); hold on;
semilogy(ListSigma,Error_2,'.-g');
semilogy(ListSigma,Error_3,'.-b');
hold off
xline(2*pi);
legend('f1','f2','f3','$2\pi$','Interpreter','latex');
xlabel("$\sigma$","Interpreter","latex");

%% Dependence of the accuracy on the epsilon
% uses the same error grid as above

ListTol = logspace(-15,-6,100);
Error_1 = []; Error_2 = []; Error_3 = [];
N_C1 = []; N_C2 = []; N_C3 = [];
for tol = ListTol
    [g, ~, C] = square_approximation(f1, [0,1,0,1], 'pole_x', [0,1], ...
        'pole_y', [0,1], 'tol', tol, 'N1', 40);
    G = g(X1,Y1);
    E = max(abs(F1-G),[],"all");
    Error_1 = [Error_1, E];
    N_C1 = [N_C1, norm(C(:))];

    [g, ~, C] = square_approximation(f2, [0,1,0,1], 'pole_x', [0], ...
        'pole_y', [0], 'tol', tol, 'N1', 40);
    G = g(X2,Y2);
    E = max(abs(F2-G),[],"all");
    Error_2 = [Error_2, E];
    N_C2 = [N_C2, norm(C,"fro")];

    [g, ~, C] = square_approximation(f3, [0,1,-1,1], 'trig_y', ...
        'pole_x', [3/4,1], 'tol', tol, 'N1', 80);
    G = g(r,t');
    E = max(abs(F3-G),[],"all");
    Error_3 = [Error_3, E];
    N_C3 = [N_C3, norm(C(:))];

    sprintf("tol = %.2e,", tol)
end

figure;
loglog(ListTol,Error_1,'.-r'); hold on;
loglog(ListTol,Error_2,'.-g'); 
loglog(ListTol,Error_3,'.-b'); xline(1e-10); 
legend('f1','f2','f3','Interpreter','latex');
ylabel('error'); xlabel('$\epsilon$','Interpreter','latex');
xlim([1e-15, 1e-6]);

%% More in-depth analysis of the influence of epsilon on the overall 
% behaviour --> comparison between epsilon = 1e-10 and 1e-12

err1_10 = []; err1_14 = []; cnorm1_10 = []; cnorm1_14 = [];
err3_10 = []; err3_14 = []; cnorm3_10 = []; cnorm3_14 = [];

Nlist = 5:10:200;

for N = Nlist
    % f1, tol = 1e-10
    [g, ~, C] = square_approximation(f1, [0,1,0,1], 'pole_x', [0,1], ...
        'pole_y', [0,1], 'N1', N, 'tol', 1e-10);
    G = g(X1,Y1);
    E = max(abs(F1-G),[],"all");
    err1_10 = [err1_10, E];
    cnorm1_10 = [cnorm1_10 norm(C(:))];
    % f1, tol = 1e-14
    [g, ~, C] = square_approximation(f1, [0,1,0,1], 'pole_x', [0,1], ...
        'pole_y', [0,1], 'N1', N, 'tol', 1e-14);
    G = g(X1,Y1);
    E = max(abs(F1-G),[],"all");
    err1_14 = [err1_14, E];
    cnorm1_14 = [cnorm1_14 norm(C(:))];

    % f3, tol = 1e-10
    [g, ~, C] = square_approximation(f3, [0,1,-1,1], 'trig_y', ...
        'pole_x', [3/4,1], 'N1', N, 'tol', 1e-10);
    G = g(r,t');
    E = max(abs(F3-G),[],"all");
    err3_10 = [err3_10, E];
    cnorm3_10 = [cnorm3_10 norm(C(:))];
    % f3, tol = 1e-14
    [g, ~, C] = square_approximation(f3, [0,1,-1,1], 'trig_y', ...
        'pole_x', [3/4,1], 'N1', N, 'tol', 1e-14);
    G = g(r,t');
    E = max(abs(F3-G),[],"all");
    err3_14 = [err3_14, E];
    cnorm3_14 = [cnorm3_14 norm(C(:))];
end

figure;
subplot(1,2,1);
semilogy(sqrt(Nlist), err1_10, '.-r'); hold on;
semilogy(sqrt(Nlist), err1_14, '--r','HandleVisibility','off');
semilogy(sqrt(Nlist), err3_10, '.-g');
semilogy(sqrt(Nlist), err3_14, '--g','HandleVisibility','off');
legend('f1','f3');
xlabel('sqrt(N1)'); ylabel('error');
axis square

subplot(1,2,2);
semilogy(sqrt(Nlist), cnorm1_10, '.-r'); hold on;
semilogy(sqrt(Nlist), cnorm1_14, '--r','HandleVisibility','off');
semilogy(sqrt(Nlist), cnorm3_10, '.-g'); 
semilogy(sqrt(Nlist), cnorm3_14, '--g','HandleVisibility','off');
legend('f1','f3');
xlabel('sqrt(N1)'); ylabel('coefficient norm');
axis square