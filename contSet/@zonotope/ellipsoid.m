function E = ellipsoid(Z,varargin)
% ellipsoid - converts a zonotope to an ellipsoid
%
% Syntax:
%    E = ellipsoid(Z)
%    E = ellipsoid(Z,mode)
%
% Inputs:
%    Z - zonotope object
%    mode - (optional) specifies whether function uses a bound on the 
%               respective zonotope norm or the exact value:
%               - 'outer:exact':    Uses priv_MVEE(Z)
%               - 'outer:norm':     Uses priv_encEllipsoid with exact norm
%                                   value
%               - 'outer:norm_bnd': Uses priv_encEllipsoid(E) with upper
%                                   bound for norm value (default)
%               - 'inner:exact':    Uses priv_MVIE(Z)
%               - 'inner:norm'      Uses priv_inscEllipsoid(E,'exact') with
%                                   exact norm value
%               - 'inner:norm_bnd': Not implemented yet, throws error
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope([1;-1],[2 -4 3 2 1; 3 2 -4 -2 1]);
%    E = ellipsoid(Z);
%
%    figure; hold on;
%    plot(Z);
%    plot(E,[1,2],'r');
%
% References:
%    [1] V. Gaßmann, M. Althoff. "Scalable Zonotope-Ellipsoid Conversions
%        using the Euclidean Zonotope Norm", 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       11-October-2019
% Last update:   05-June-2021 (MA, include degenerate case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,2);

% set default method
mode = setDefaultValues({'outer:norm_bnd'},varargin);

% obtain information about zonotope
n = dim(Z);

% zonotope is just a point
if representsa_(Z,'point',eps)
    E = ellipsoid(zeros(n),Z.c);
    return
end

% exact information
c = Z.c; G = Z.G;
Grank = rank(G);

% reduce dimension of zonotope if degenerate
if n ~= Grank
    % compute QR decomposition
    [Q, R] = qr(G);
    
    % obtain reduced zonotope centered around the origin
    Z = zonotope([zeros(Grank,1), R(1:Grank,:)]);
    
    % obtain transformation matrix
    T = Q(:,1:Grank);
end

% Zonotope is a parallelotope 
if size(Z.G,2) == n
    fac = n;
    if startsWith(mode,'inner')
        fac = 1;
    end
    E = ellipsoid(fac*(G*G'),Z.c);
% Zonotope is not a parallelotope
else
    switch mode
        % Exact outer enclosure
        case 'outer:exact'
            E = priv_MVEE(Z);
        % Norm-based outer enclosure
        case 'outer:norm'
            E = priv_encEllipsoid(Z,'exact');
        % Norm-based outer enclosure with upper bound for norm value
        case 'outer:norm_bnd'
            E = priv_encEllipsoid(Z,'ub_convex');
        % Exact inner enclosure
        case 'inner:exact'
            E = priv_MVIE(Z);
        % Norm-based inner enclosure
        case 'inner:norm'
            E = priv_inscEllipsoid(Z);
        %case 'inner:norm_bnd' % we do not implement the test used to compute the
        %necessary lower bound on Z as in [1] since this test generally does
        %not result in a very good lower bound
        %    E = priv_inscEllipsoid(Z);
        otherwise
             throw(CORAerror('CORA:wrongValue','second',...
                 "'outer:exact','outer:norm','outer:norm_bnd','inner:exact' or 'inner:norm'"));
    end
end

% convert to original space if zonotope is not full-dimensional
if n ~= Grank
    E = T*E + c;
end

% ------------------------------ END OF CODE ------------------------------
