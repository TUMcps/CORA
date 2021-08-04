function Z = zonotope(E,m,mode)
% zonotope - overapproximates an ellipsoid by a zonotope
%
% Syntax:  
%    E = enc_zonotope(Z,m,mode)
%
% Inputs:
%    E       - ellipsoid object
%    m       - number of generators
%    mode    - (Optional) Specifies whether function uses a lower bound on the 
%               minimum zonotope norm or the exact value:
%               * 'o:box':     overapprox. parallelotope using
%                              enc_parallelotope
%               * 'o:norm':    Uses enc_zonotope(E,m,'exact') with exact norm value
%               * 'o:norm:bnd':Not implemented yet, throws error
%               * 'i:box':     inner approx. parallelotope using
%                              insc_parallelotope
%               * 'i:norm'     Uses insc_zonotope(E,m,'exact') with exact norm
%                              value
%               * 'i:norm:bnd':Uses insc_zonotope(E,m) with an bound on the
%                              norm value
%               * default:     same as 'o:norm:bnd'
%
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid.generateRandom(0,2);
%    Zenc = zonotope(E,10,'o:norm');
%    Zinsc = zonotope(E,10,'i:norm');
%    Zbox = zonotope(E);
%    plot(E);
%    hold on
%    plot(Zinsc);
%    plot(Zenc);
%    plot(Zbox);
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

% Author:       Victor Gassmann
% Written:      11-October-2019
% Last update:  08-June-2021 (handle degenerate case here, remove from
%                             sub-files)
% Last revision:---

%------------- BEGIN CODE --------------
default = 'o:norm:bnd';
if ~exist('m','var')
    if ~exist('mode','var')
        mode = 'o:box';
    else
        mode = default;
    end
end

% compute rank and dimension of ellipsoid
rankE = rank(E);
dimE = dim(E);
c = center(E);

% handle degenerate case
if ~(rankE == dimE)
    if rankE==0
        Z = zonotope(E.q);
        return;
    end
    % compute svd to find U matrix transforming E.Q to diagonal matrix (to
    % isolate degenerate dimensions)
    [U,Qt,~] = svd(E.Q);
    
    % construct non-degenerate ellipsoid
    E = ellipsoid(Qt(1:rankE,1:rankE));
    
    % construct revert transformation matrix
    T = U(:,1:rankE);
end


switch mode
    case 'o:box'
        Z = enc_parallelotope(E);
    case 'o:norm'
        Z = enc_zonotope(E,m,'exact');
    %case 'o:norm:bnd'% we do not implement the test used to compute the
    %necessary lower bound on Z as in [1] since this test generally does
    %not result in a very good lower bound
    %    Z = enc_zonotope(E,m);
    case 'i:box'
        Z = insc_parallelotope(E);
    case 'i:norm'
        Z = insc_zonotope(E,m,'exact');
    case 'i:norm:bnd'
        Z = insc_zonotope(E,m);
    otherwise
        error('Wrong value for argument "mode"');
end

if ~(rankE == dimE)
    Z = T*Z + c;
end
%------------- END OF CODE --------------