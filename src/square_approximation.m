function varargout = square_approximation(f, varargin)
% g = square_approximation(F, [A B C D]) returns a multivariate rational function of
% the form g = Phi(x)*C*Psi(y) on the rectangle [-1,1] x [-1,1].
%
% [g, res] = square_approximation(F, [A B C D]) returns a multivariate rational function 
% on the rectangle [A B] x [C x D], along with least-square residual.
%
% [g, res, C] square_approximation(F, 'trig') constructs a multivariate rational
% approximation with a bivariate Fourier polynomial part with least-square residual
% and coefficient matrix C.
%
% g = square_approximation(F, 'trigx') constructs a Fourier x Chebyshev expansion 
% that is periodic in x only.
%
% g = square_approximation(F, 'trigy') constructs a Chebyshev x Fourier expansion 
% that is periodic in y only. 

% Parse the input
[dom, trig_x, trig_y, pole_x, pole_y, N1, N2, N3, M, sigma, tol] = ParseInputs(varargin{:});

% Define domain
dom_x = dom(1:2);
dom_y = dom(3:4);

% Construct grid to evaluate the function
X = [];
if not(isempty(pole_x))
    X_cluster = logspace(log10(eps), 0, M)';
    for p = pole_x 
        X = [X; (dom_x(1)-p)*X_cluster+p; (dom_x(2)-p)*X_cluster+p];
    end
end

% Add trigonometric or Chebyshev points
if trig_x
    X = [X; trigpts(2*N2, dom_x)];
else
    X = [X; chebpts(2*N2, dom_x)];
end

% Remove duplicates
X = unique(X);

% Do the same for y
Y = [];
if not(isempty(pole_y))
    Y_cluster = logspace(log10(eps), 0, M)';
    for p = pole_y
        Y = [Y; (dom_y(1)-p)*Y_cluster+p; (dom_y(2)-p)*Y_cluster+p];
    end
end
if trig_y
    Y = [Y; trigpts(2*N3, dom_y)];
else
    Y = [Y; chebpts(2*N3, dom_y)];
end

% Remove duplicates
Y = unique(Y);

% Evaluate the function at the grid of points
[x,y] = ndgrid(X,Y);
F = f(x(:),y(:));
F = reshape(F, [length(X), length(Y)]);

% Define poles
pj = -exp(-sigma*(sqrt(N1)-sqrt(1:N1)));

% Rational + Chebyshev polynomial part on x

% Define poles x
if not(isempty(pole_x))
    pj_x = [];

    % First order poles
    for p = pole_x
        pj_x = [pj_x, p*(1+eps)+1i*pj, p*(1+eps)-1i*pj];
    end
    pj_x_repeat = repmat(pj,1,2*length(pole_x));
    r_x = @(x) pj_x_repeat ./ (x-pj_x);
else
    r_x = @(x) [];
end

if trig_x
    px = trigpoly(-N2:N2, dom_x);
else
    px = chebpoly(0:N2, dom_x);
end
Phi = @(x)[r_x(x), px(x)];
Phi_x = Phi(X);
% Rational + Chebyshev polynomial part on y

% Define poles y
if not(isempty(pole_y))
    pj_y = [];
    for p = pole_y
        pj_y = [pj_y, p*(1+eps)+1i*pj, p*(1+eps)-1i*pj];
    end
    pj_y_repeat = repmat(pj,1,2*length(pole_y));
    r_y = @(x) pj_y_repeat ./ (x-pj_y);
else
    r_y = @(x) [];
end

if trig_y
    py = trigpoly(-N3:N3, dom_y);
else
    py = chebpoly(0:N3, dom_y);
end
Psi = @(x)[r_y(x), py(x)];
Psi_y = Psi(Y);

% Solve linear system Phi*C*Psi' = F
C = truncated2_svd_reg(Phi_x, Psi_y, F, tol);
% Compute linear solver error:
res = max(abs(Phi_x*C*Psi_y.'-F), [], "all");
sprintf("Least-square error = %.2e", res)

% Return approximant
g = @(x,y) Phi(x)*C*Psi(y).';

varargout{1} = g;
if nargout >= 2
    varargout{2} = res;
end
if nargout == 3
    varargout{3} = C;
end
end

function [dom, trig_x, trig_y, pole_x, pole_y, N1, N2, N3, M, sigma, tol] = ParseInputs(varargin)

% Initialize value of parameters
dom = [-1,1,-1,1];
trig_x = false;
trig_y = false;
pole_x = [];
pole_y = [];

% Degree of rational expansion
N1 = 150;

% Degree of smooth part in the x coordinates
N2 = ceil(1.3*sqrt(N1));

% Degree of smooth part in the y coordinates
N3 = ceil(1.3*sqrt(N1));

% Tolerance for solving least-square problem
tol = 1e-14;

% Pole spacing
sigma = 2*pi;

% Check if number of sample points is provided
is_M = false;
args = varargin;

% Try to parse out the domain
if ( ~isempty(args) )
    if ( isnumeric(args{1}) && (length(args{1}) == 4) )
        dom = args{1};
        args(1) = [];
    end
end

% Parse out other parameters
while ( ~isempty(args) )

    if strcmpi(args{1}, 'trig_x')
        trig_x = true;
        args(1) = [];
    elseif strcmpi(args{1}, 'trig_y')
        trig_y = true;
        args(1) = [];
    elseif strcmpi(args{1}, 'trig')
        trig_x = true;
        trig_y = true;
        args(1) = [];
    elseif strcmpi(args{1}, 'pole_x')
        if length(args) > 1 && isnumeric(args{2})
            pole_x = args{2};
            args(1:2) = [];
        else
            error('Could not parse pole x locations.');
        end
    elseif strcmpi(args{1}, 'pole_y')
        if length(args) > 1 && isnumeric(args{2})
            pole_y = args{2};
            args(1:2) = [];
        else
            error('Could not parse pole y locations.');
        end
    elseif strcmpi(args{1}, 'N1')
        if length(args) > 1 && isnumeric(args{2})
            N1 = args{2};
            args(1:2) = [];
        else
            error('Could not parse N1.');
        end
    elseif strcmpi(args{1}, 'N2')
        if length(args) > 1 && isnumeric(args{2})
            N2 = args{2};
            args(1:2) = [];
        else
            error('Could not parse N2.');
        end
    elseif strcmpi(args{1}, 'N3')
        if length(args) > 1 && isnumeric(args{2})
            N3 = args{2};
            args(1:2) = [];
        else
            error('Could not parse N3.');
        end
    elseif strcmpi(args{1}, 'M')
        if length(args) > 1 && isnumeric(args{2})
            M = args{2};
            is_M = true;
            args(1:2) = [];
        else
            error('Could not parse M.');
        end
    elseif strcmpi(args{1}, 'sigma')
        if length(args) > 1 && isnumeric(args{2})
            sigma = args{2};
            args(1:2) = [];
        else
            error('Could not parse sigma.');
        end
    elseif strcmpi(args{1}, 'tol')
        if length(args) > 1 && isnumeric(args{2})
            tol = args{2};
            args(1:2) = [];
        else
            error('Could not parse tol.');
        end
    else
        error('Could not parse input argument sequence.');
    end
end

% Set number of sample points
if not(is_M)
    M = 5*N1;
end

end
