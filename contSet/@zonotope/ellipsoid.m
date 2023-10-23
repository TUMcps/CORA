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
%               - 'outer:exact':    Uses MVEE(Z)
%               - 'outer:norm':     Uses enc_ellipsoid(E,'exact') with
%                                   exact norm value
%               - 'outer:norm_bnd': Uses enc_ellipsoid(E) with upper bound
%                                   for norm value (default)
%               - 'inner:exact':    Uses MVIE(Z)
%               - 'inner:norm'      Uses insc_ellipsoid(E,'exact') with
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
%           using the Euclidean Zonotope Norm", 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       11-October-2019
% Last update:   05-June-2021 (MA, include degenerate case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default method
mode = setDefaultValues({'outer:norm_bnd'},varargin);

% obtain rank of zonotope
n = dim(Z);
G = Z.G;
c = Z.c;
Grank = rank(G);


% zonotope is the origin
if representsa_(Z,'origin',eps)

    E = ellipsoid(zeros(n),zeros(n,1));
    return

% reduce dimension of zonotope if not full dimensional
elseif n ~= Grank
    
    % compute QR decomposition
    [Q, R] = qr(G);
    
    % obtain reduced zonotope centered around the origin
    Z = zonotope([zeros(Grank,1), R(1:Grank,:)]);
    
    % obtain transformation matrix
    T = Q(:,1:Grank);
end

if size(Z.G,2) == n
    fac = n;
    if startsWith(mode,'inner')
        fac = 1;
    end
    E = ellipsoid(fac*(G*G'),Z.c);
else
    switch mode
        case 'outer:exact'
            E = MVEE(Z);
        case 'outer:norm'
            E = enc_ellipsoid(Z,'exact');
        case 'outer:norm_bnd'
            E = enc_ellipsoid(Z);
        case 'inner:exact'
            E = MVIE(Z);
        case 'inner:norm'
            E = insc_ellipsoid(Z,'exact');
        %case 'inner:norm_bnd' % we do not implement the test used to compute the
        %necessary lower bound on Z as in [1] since this test generally does
        %not result in a very good lower bound
        %    E = insc_ellipsoid(Z);
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
