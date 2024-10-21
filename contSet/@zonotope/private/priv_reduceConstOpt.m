function Zred = priv_reduceConstOpt(Z,order,method,alg)
% priv_reduceConstOpt - method to reduce the order of a zonotope by 
%    constraint convex optimization: minimize volume of the zonotope used
%    for overapproximation abs(det(C)) subject to all points of the
%    zonotope are inside sum(abs(C^-1 * G)) <= 1
%
% Syntax:
%    Zred = priv_reduceConstOpt(Z,order,method,alg)
%
% Inputs:
%    Z - zonotope object 
%    order - order of the reduced zonotope (order=#generators/dim)
%    method - minimize det or Frobenius norm
%    alg - algorithm used by the solver fmincon
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Authors:       Anna Kopetzki, Matthias Althoff
% Written:       11-September-2016 (AK)
% Last update:   27-June-2018 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize reduced zonotope
Zred = Z;

% pick generators to reduce
[~, Gunred, Gred] = pickedGenerators(Z,order);

if isempty(Gred)
    return
end

% Initialize the transformation matrix C0 (overapproximation using PCA)
Zred = priv_reducePCA(Z,1);
C0 = generators(Zred);
%C0=pcaInit(G);
%     C = C0;

% Nonlinear constraints, fmincon with more iterations 
options = optimoptions(@fmincon,'Algorithm', alg,...
    'MaxIterations',5000,...
    'MaxFunctionEvaluations',100000,...
    'Display','off');

switch method
    case 'det'
        Gred = fmincon(@(X)minVolEst(X,Z.G),C0,[],[],[],[],[],[], @(X)zonoConst(X,Z.G),options);
    case 'frob'
        Gred = fmincon(@(X)minVolApprox(X,Z.G),C0,[],[],[],[],[],[], @(X)zonoConst(X,Z.G),options);
    case 'qr'
        Gred = aux_qr(C0,Z.G,options);
    case 'svd'
        Gred = aux_svd(C0,Z.G,options);
end

% build reduced zonotope
Zred.c = Z.c;
Zred.G = [Gunred,Gred];

end


% Auxiliary functions -----------------------------------------------------

function Gred = aux_qr(C0,G,options)

n = size(G,1);
[Q0, R0] = qr(C0);
X0 = [reshape(Q0, n*n, 1); reshape(R0, n*n, 1)];
Y_vec = fmincon(@(X)qrLogVol(X,G),X0,[],[],[],[],[],[], @(X)zonoQRConst(X,G),options);
Y = reshape(Y_vec, n, 2*n);
Q = Y(:,1:n);
R = triu(Y(:,(n+1):2*n));
Gred = Q*triu(R);

end

function Gred = aux_svd(C0,G,options)

n = size(G,1);
[U0, S0, V0] = svd(C0);
X0 = [reshape(U0, n*n, 1); reshape(S0, n*n, 1); reshape(V0, n*n, 1)];
Y_vec = fmincon(@(X)aux_svdLogVol(X,G),X0,[],[],[],[],[],[], @(X)aux_zonoSVDConst(X,G),options);
Y = reshape(Y_vec, n, 3*n);
U = Y(:,1:n);
S = Y(:,(n+1):2*n);
V = Y(:,(2*n+1):3*n);
Gred = U * diag(diag(S)) * V';

end

function vol = aux_svdLogVol(X, G)
    n = size(G,1);
    Y = reshape(X, n, 3*n);
    S = Y(:,(n+1):2*n);
    
    vol = sum(log(diag(abs(S)))); % log(diag(abs(S))) can be < 0
end

function [c, ceq] = aux_zonoSVDConst(X, G)
    n = size(G,1);
    Y = reshape(X, n, 3*n);
    U = Y(:,1:n);
    S = Y(:,(n+1):2*n);
    V = Y(:,(2*n+1):3*n);
    
    C_inv = V * diag(diag(1 ./S)) * U'; 
    c = sum(abs(C_inv * G), 2) - ones(n,1);
    ceq = [U*U' - diag(ones(1,n)); V*V' - diag(ones(1,n))];

end

% ------------------------------ END OF CODE ------------------------------
