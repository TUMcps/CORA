function [E] = ellipsoid(Z,mode)
% enc_ellipsoid - Overapproximates a zonotope by an ellipsoid
%
% Syntax:  
%    E = ellipsoid(Z,comptype)
%
% Inputs:
%    Z - zonotope object
%    mode - (optional) specifies whether function uses a bound on the 
%               respective zonotope norm or the exact value:
%               - 'o:exact':   Uses MVEE(Z)
%               - 'o:norm':    Uses enc_ellipsoid(E,'exact') with exact norm value
%               - 'o:norm:bnd':Uses enc_ellipsoid(E) with upper bound for norm
%                              value
%               - 'i:exact':   Uses MVIE(Z)
%               - 'i:norm'     Uses insc_ellipsoid(E,'exact') with exact norm
%                              value
%               - 'i:norm:bnd':Not implemented yet, throws error
%               - default:     same as 'o:norm:bnd'
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope(rand(2,5));
%    E = ellipsoid(Z);%same as ellipsoid(Z,'o:norm:bnd')
%    plot(Z);
%    hold on
%    plot(E);
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

% Author:       Victor Gassmann, Matthias Althoff
% Written:      11-October-2019
% Last update:  05-June-2021 (MA, include degenerate case)
% Last revision:---

%------------- BEGIN CODE --------------
default = 'o:norm:bnd';
if ~exist('mode','var')
    mode = default;
end

% obtain rank of zonotope
Zdim = dim(Z);
G = generators(Z);
c = center(Z);
Grank = rank(G);


% reduce dimension of zonotope if not full dimensional
if ~(Zdim == Grank)
    
    % compute QR decomposition
    [Q, R] = qr(G);
    
    % obtain reduced zonotope centered around the origin
    Z = zonotope([zeros(Grank,1), R(1:Grank,:)]);
    
    % obtain transformation matrix
    T = Q(:,1:Grank);
end

if size(generators(Z),2)==dim(Z)
    n = dim(Z);
    G = generators(Z);
    fac = n;
    if startsWith(mode,'i')
        fac = 1;
    end
    E = ellipsoid(fac*(G*G'),center(Z));
else
    switch mode
        case 'o:exact'
            E = MVEE(Z);
        case 'o:norm'
            E = enc_ellipsoid(Z,'exact');
        case 'o:norm:bnd'
            E = enc_ellipsoid(Z);
        case 'i:exact'
            E = MVIE(Z);
        case 'i:norm'
            E = insc_ellipsoid(Z,'exact');
        %case 'i:norm:bnd' % we do not implement the test used to compute the
        %necessary lower bound on Z as in [1] since this test generally does
        %not result in a very good lower bound
        %    E = insc_ellipsoid(Z);
        otherwise
            error('Wrong value for argument "mode"');
    end
end

% convert to original space if zonotope is not full dimensional
if ~(Zdim == Grank)
    E = T*E + c;
end
%------------- END OF CODE --------------