%% Less stable alternative for the 2D truncated SVD solver
% See Algorithm A.1 in the paper

function x = truncated2_svd(A, B, F, rtol)
    [U1,S1,V1] = svd(A,"econ");
    singvals1 = diag(S1);
    K = find(singvals1/singvals1(1) > rtol, 1, 'last');
    S1_eps_inv = zeros(size(S1));
    S1_eps_inv(1:K,1:K) = diag(singvals1(1:K).^(-1));
    
    [U2,S2,V2] = svd(B.',"econ");
    singvals2 = diag(S2);
    K = find(singvals2/singvals2(1) > rtol, 1, 'last');
    S2_eps_inv = zeros(fliplr(size(S2)));
    S2_eps_inv(1:K,1:K) = diag(singvals2(1:K).^(-1));
    x = V1*(S1_eps_inv'*(U1'*F*V2)*S2_eps_inv)*U2';
end
