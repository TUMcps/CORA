function Z = zonotope(E,varargin)
% zonotope - converts an ellipsoid to a zonotope
%
% Syntax:
%    Z = zonotope(E)
%    Z = zonotope(E,mode)
%    Z = zonotope(E,mode,nrGen)
%
% Inputs:
%    E - ellipsoid object
%    mode - (optional) Specifies whether function uses a lower bound on the 
%           minimum zonotope norm or the exact value:
%           * 'outer:box':      overapprox. parallelotope using
%                               priv_encParallelotope
%           * 'outer:norm':     uses priv_encZonotope with exact norm value
%           * 'outer:norm_bnd': not implemented yet (throws error)
%           * 'inner:box':      inner approx. parallelotope using
%                               priv_inscParallelotope
%           * 'inner:norm'      uses priv_inscZonotope with exact norm value
%           * 'inner:norm_bnd': uses priv_inscZonotope with an bound on the
%                               norm value
%           * default:          same as 'outer:norm_bnd'
%    nrGen - (optional) number of generators
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z_enc = zonotope(E,'outer:norm',10);
%    Z_insc = zonotope(E,'inner:norm',10);
%    Z_box = zonotope(E);
%
%    figure; hold on;
%    plot(E);
%    plot(Z_insc,[1,2],'r');
%    plot(Z_enc,[1,2],'k');
%    plot(Z_box,[1,2],'m');
%
% References:
%    [1] V. Gaßmann, M. Althoff. "Scalable Zonotope-Ellipsoid Conversions
%        using the Euclidean Zonotope Norm", 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_encZonotope, priv_encParallelotope, priv_inscZonotope,
%    priv_inscParallelotope

% Authors:       Victor Gassmann
% Written:       11-October-2019
% Last update:   08-June-2021 (moved handling of degenerate case here)
%                04-July-2022 (VG, class array cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,3);

[mode,nrGen] = setDefaultValues({'outer:box',dim(E)},varargin);

% before v2025, one could call the conversion with only the number of
% generators, which has been removed to simplify the interface
if isnumeric(mode)
    CORAwarning("CORA:interface","ellipsoid/zonotope","CORA v2025.0.0");
end

inputArgsCheck({{E,'att','ellipsoid'};
                {mode,'str',{'outer:box','outer:norm','outer:norm_bnd',...
                             'inner:box','inner:norm','inner:norm_bnd'}}; ...
                {nrGen,'att','numeric',{'scalar','positive','integer'}}});
% nrGen not used if outer:box or inner:box since there the number of
% generators is fixed!

% handle empty case
if representsa(E, 'emptySet')
    Z = zonotope.empty(dim(E));
    return
end

% compute rank and dimension of ellipsoid
rankE = rank(E);
c = center(E);
isDeg = rankE ~= dim(E);

% handle degenerate case
if isDeg
    % ellipsoid is just a point
    if rankE == 0
        Z = zonotope(c);
        return;
    end

    % compute svd to find U matrix transforming the shape matrix to a 
    % diagonal matrix (to isolate degenerate dimensions)
    [U,Qt,~] = svd(E.Q);
    
    % construct non-degenerate ellipsoid
    E = ellipsoid(Qt(1:rankE,1:rankE));
    
    % construct revert transformation matrix
    T = U(:,1:rankE);
end


switch mode
    case 'outer:box'
        Z = priv_encParallelotope(E);
    case 'outer:norm'
        Z = priv_encZonotope(E,nrGen);
    case 'outer:norm:bnd'
        throw(CORAerror('CORA:notSupported',"mode = 'outer:norm:bnd'"));
        % we do not implement the test used to compute the
        % necessary lower bound on Z as in [1] since this test generally does
        % not result in a very good lower bound
    %    Z = priv_encZonotope(E,nrGen,'approx');
    case 'inner:box'
        Z = priv_inscParallelotope(E);
    case 'inner:norm'
        Z = priv_inscZonotope(E,nrGen,'exact');
    case 'inner:norm_bnd'
        Z = priv_inscZonotope(E,nrGen,'ub_convex');
end

% in degenerate case, lift lower-dimensional non-degenerate ellipsoid
if isDeg
    Z = T*Z + c;
end

% ------------------------------ END OF CODE ------------------------------
