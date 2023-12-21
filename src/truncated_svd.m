%% 1D Truncated SVD solver

function x = truncated_svd(A, b, rtol)
    [U,S,V] = svd(A, 0);
    singvals = diag(S);
    K = find(singvals/singvals(1) > rtol, 1, 'last');
    S_eps_inv = zeros(size(S))';
    S_eps_inv(1:K,1:K) = diag(singvals(1:K).^(-1));
    x = V*(S_eps_inv*(U'*b));
end
