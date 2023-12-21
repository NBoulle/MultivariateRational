%function C = linearize_tensorproduct(A, B)
%
%   Make a (dense) matrix C that represents the tensor product A x B,
%   of the form [A*B(:,1) A*B(:,2) ...].
function C = linearize_tensorproduct(A, B)
    [m,n] = size(A);
    [p,q] = size(B);

    assert(m==p, "incompatible dimensions")
    C = zeros(m, n*q);
    for i = 1:q
        C(:,(i-1)*n+1:i*n) = A .* repmat(B(:,i), 1, n);
    end
end
