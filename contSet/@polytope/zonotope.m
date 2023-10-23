function Z = zonotope(P,varargin)
% zonotope - Overapproximates a polytope by a zonotope
%
% Syntax:
%    E = zonotope(P)
%    E = zonotope(P,mode)
%
% Inputs:
%    P - polytope object
%    mode - type of approximation
%              'outer': outer-approximation
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    P = polytope(rand(2,5));
%    E = zonotope(P);
%
%    figure; hold on
%    plot(P);
%    plot(E,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       17-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

mode = setDefaultValues({'outer'},varargin);

if strcmp(mode,'outer')
    % compute maximum-volume ellipsoid to get "shape" of P and center
    n = dim(P);
    E = ellipsoid(P,'inner');
    c = center(E);
    Q = E.Q;
    Q_r = sqrtm(Q);
    
    % shift and transform polytope so that it is roughly a hyperbox
    Mt = P;
    Mt.A = P.A*Q_r;
    Mt.b = P.b -P.A*c;
    if ~isempty(P.Ae)
        Mt.Ae = P.Ae*Q_r;
        Mt.be = P.be -M.Ae*c;
    end

    % compute enclosing hyperbox
    Supp_u = zeros(n,1);
    Supp_l = zeros(n,1);
    I = eye(n);
    for i=1:n
        Supp_u(i) = supportFunc(Mt,I(:,i),'upper');
        Supp_l(i) = supportFunc(Mt,I(:,i),'lower');
    end
    ct = 1/2*(Supp_u+Supp_l);
    rt = 1/2*(Supp_u-Supp_l);

    % form result
    Z = Q_r*zonotope(ct,diag(rt)) +c;
else
    throw(CORAerror("CORA:notSupported",'mode="inner" is not supported.'));
end

% ------------------------------ END OF CODE ------------------------------
