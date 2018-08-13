function [Zred]=reduceConstOpt(Z,order, method, alg)
% reduceConstOpt - method to reduce the order of a zonotope
% Method: Constraint Convex Optimization 
% Minimize volume of the zonotope used for overapproximation abs(det(C))
% Subect to all Points of the Zonotope are inside 
% sum(abs(C^-1 * G)) <= 1
%
% Syntax:  
%    [Zred]=reduceConstOpt(Z,order, method, alg)
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
% See also: 

% Author:       Anna Kopetzki, Matthias Althoff
% Written:      11-September-2016 (AK)
% Last update:  27-June-2018 (MA)
% Last revision:---

%------------- BEGIN CODE --------------

% initialize Z_red
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGenerators(Z,order);

if ~isempty(Gred)

    % get Z-matrix, center cen, generator matrix G, dimension dim from zonotope
    % Z, number of generators nrG, desNrGen = desired number of generators
    Zmatrix=get(Z,'Z');
    cen=Zmatrix(:,1);
    G=Zmatrix(:,2:end);
    dim=length(cen);

    % Initialize the transformation matrix C0 (overapproximation using PCA)
    Zred = reducePCA(Z,1);
    C0 = generators(Zred);
    %C0=pcaInit(G);
    C = C0;

    % Nonlinear constraints, fmincon with more iterations 
    options = optimoptions(@fmincon,'Algorithm', alg, 'MaxIterations',5000, 'MaxFunctionEvaluations', 100000, 'Display', 'off');


    if strcmp(method, 'det')
        Gred = fmincon(@(X)minVolEst(X,G),C0,[],[],[],[],[],[], @(X)zonoConst(X,G),options);
    elseif strcmp(method, 'frob')
        Gred = fmincon(@(X)minVolApprox(X,G),C0,[],[],[],[],[],[], @(X)zonoConst(X,G),options);
    elseif strcmp(method, 'qr')
        [Q0, R0] = qr(C0);
        X0 = [reshape(Q0, dim*dim, 1); reshape(R0, dim*dim, 1)];
        Y_vec = fmincon(@(X)qrLogVol(X,G),X0,[],[],[],[],[],[], @(X)zonoQRConst(X,G),options);
        Y = reshape(Y_vec, dim, 2*dim);
        Q = Y(:,1:dim);
        R = triu(Y(:,(dim+1):2*dim));
        Gred = Q*triu(R);
    elseif strcmp(method, 'svd')
        [U0, S0, V0] = svd(C0);
        X0 = [reshape(U0, dim*dim, 1); reshape(S0, dim*dim, 1); reshape(V0, dim*dim, 1)];
        Y_vec = fmincon(@(X)svdLogVol(X,G),X0,[],[],[],[],[],[], @(X)zonoSVDConst(X,G),options);
        Y = reshape(Y_vec, dim, 3*dim);
        U = Y(:,1:dim);
        S = Y(:,(dim+1):2*dim);
        V = Y(:,(2*dim+1):3*dim);
        Gred = U * diag(diag(S)) * V';
    end
end

%build reduced zonotope
Zred.Z=[center,Gunred,Gred];



function vol=svdLogVol(X, G)

dim = size(G,1);
Y = reshape(X, dim, 3*dim);
S = Y(:,(dim+1):2*dim);

vol = sum(log(diag(abs(S)))); % log(diag(abs(S))) can be < 0


function [c, ceq]=zonoSVDConst(X, G)

dim = size(G,1);
Y = reshape(X, dim, 3*dim);
U = Y(:,1:dim);
S = Y(:,(dim+1):2*dim);
V = Y(:,(2*dim+1):3*dim);


C_inv = V * diag(diag(1 ./S)) * U'; 
c = sum(abs(C_inv * G), 2) - ones(dim,1);
ceq = [U*U' - diag(ones(1,dim)); V*V' - diag(ones(1,dim))];

%-------------------- END CODE --------------------
