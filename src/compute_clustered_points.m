function [X,Y] = compute_clustered_points(f_curve, N_tangent, N_normal, dom, threshold)
% [X,Y] = compute_clustered_points(f_curve, N_tangent, N_normal, dom[, threshold])
% 
% returns coordinates of points clustering to the curve parameterized by f_curve.
% f_curve is a 2D function handle vanishing at the singularity curve
% N_tangent is the number of points along the curve
% N_normal is the number of clustering points in the normal direction
% dom = [xmin, xmax, ymin, ymax] is the domain size
% Optionally, threshold specifies how close the points can be to the curve.

if nargin == 4
    threshold = eps;
end

% Compute the zero level curve
f = chebfun2(f_curve, dom);
g = roots(f);

% Compute gradient of f
grad_f = grad(f);

% Discretization along the tangential direction
t_disc = linspace(-1,1,N_tangent);

% Clustering points
X_cluster = logspace(log10(threshold),0,N_normal);

% Define arrays of points
X = [];
Y = [];

% Loop over the number of disconnected components
for i = 1:size(g,2)

    % Loop over the number of points
    for j = 1:length(t_disc)

        % Get coordinate of points in the curve
        t = t_disc(j);
        xint = real(g(t,i));
        yint = imag(g(t,i));
        
        % Compute the normal at that point
        n = [real(grad_f(xint,yint)), imag(grad_f(xint,yint))];
        n = n(:,1)/norm(n);
        
        % Compute interestion between the normal equation and circle centered 
        % at (xint,yint) and radius 0.5
        t_int = 0.5/sqrt(n(1)^2+n(2)^2);
        xout = [-n(1)*t_int+xint, n(1)*t_int+xint];
        yout = [-n(2)*t_int+yint, n(2)*t_int+yint];

        % Rescale sample points
        X1 = ((xout(1)-xint)+1i*(yout(1)-yint))*X_cluster+xint+1i*yint;
        X2 = ((xout(2)-xint)+1i*(yout(2)-yint))*X_cluster+xint+1i*yint;

        % Select points inside the domain
        Z = [X1, X2];
        Z = Z(real(Z) >= dom(1));
        Z = Z(real(Z) <= dom(2));
        Z = Z(imag(Z) >= dom(3));
        Z = Z(imag(Z) <= dom(4));
        
        % Add points
        X = [X, real(Z)];
        Y = [Y, imag(Z)];
    end
end

end