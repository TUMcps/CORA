function [E] = ellipsoid(Z,mode)
% enc_ellipsoid - Overapproximates a zonotope by an ellipsoid
%
% Syntax:  
%    E = ellipsoid(Z,comptype)
%
% Inputs:
%    Z       - zonotope object
%    mode    - (Optional) Specifies whether function uses a bound on the 
%               respective zonotope norm or the exact value:
%               * 'o:exact':   Uses MVEE(Z)
%               * 'o:norm':    Uses enc_ellipsoid(E,'exact') with exact norm value
%               * 'o:norm:bnd':Uses enc_ellipsoid(E) with exact norm
%                              value
%               * 'u:exact':   Uses MVIE(Z)
%               * 'u:norm'     Uses insc_ellipsoid(E,'exact') with exact norm
%                              value
%               * 'u:norm:bnd':Not implemented yet, throws error
%               * default:     same as 'o:norm:bnd'
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
if ~exist('mode','var')
    mode = default;
end

switch mode
    case 'o:exact'
        E = MVEE(Z);
    case 'o:norm'
        E = enc_ellipsoid(Z,'exact');
    case 'o:norm:bnd'
        E = enc_ellipsoid(Z);
    case 'u:exact'
        E = MVIE(Z);
    case 'u:norm'
        E = insc_ellipsoid(Z,'exact');
    case 'u:norm:bnd'
        E = insc_ellipsoid(Z);
    otherwise
        error('Wrong value for argument "mode"');
end
%------------- END OF CODE --------------