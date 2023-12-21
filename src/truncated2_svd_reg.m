%% Solve problem AXB' = F using regularized truncated svd
% See Algorithm 2.1 in the paper

function x = truncated2_svd_reg(A, B, F, rtol)
    [U1,S1,V1] = svd(A,"econ");
    [U2,S2,V2] = svd(B,"econ");
    invSV = diag(S1).^(-1) * (diag(S2).^(-1))';
    C = U1'*(F*(U2').');
    invSV(invSV > 1/(S1(1,1)*S2(1,1)*rtol)) = 0;
    x = V1*(C.*invSV)*V2.';
end
