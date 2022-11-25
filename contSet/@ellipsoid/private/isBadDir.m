function res = isBadDir(L,E1,E2)
    [U1,S1,~] = svd(E1.Q);
    [U2,~,~] = svd(E2.Q);
    T = U2'*inv(sqrtm(S1))*U1';
    r = 1/max(diag(T*E2.Q*T'));
    res = true(1,size(L,2));
    for i=1:size(L,2)
        l = L(:,i);
        res(i) = sqrt(l'*E1.Q*l)/sqrt(l'*E2.Q*l) > r;
    end
end
