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
%               * 'u:box':     underapprox. parallelotope using
%                              insc_parallelotope
%               * 'u:norm'     Uses insc_zonotope(E,m,'exact') with exact norm
%                              value
%               * 'u:norm:bnd':Uses insc_zonotope(E,m) with an bound on the
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
%    Zinsc = zonotope(E,10,'u:norm');
%    Zbox = zonotope(E);
%    plot(E);
%    hold on
%    plot(Zinsc);
%    plot(Zenc);
%    plot(Zbox);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      11-October-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
default = 'o:norm:bnd';
if ~exist('m','var')
    if ~exist('mode','var')
        Z = enc_parallelotope(E);
        return;
    else
        mode = default;
    end
end

switch mode
    case 'o:box'
        Z = enc_parallelotope(E);
    case 'o:norm'
        Z = enc_zonotope(E,m,'exact');
    case 'o:norm:bnd'
        Z = enc_zonotope(E,m);
    case 'u:box'
        Z = insc_parallelotope(E);
    case 'u:norm'
        Z = insc_zonotope(E,m,'exact');
    case 'u:norm:bnd'
        Z = insc_zonotope(E,m);
    otherwise
        error('Wrong value for argument "mode"');
end
%------------- END OF CODE --------------