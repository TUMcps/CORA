function [res,S] = representsa_(zB,type,tol,varargin)
% representsa_ - checks if a polynomial zonotope can also be represented by
%    a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(zB,type,tol)
%    [res,S] = representsa_(zB,type,tol)
%
% Inputs:
%    zB - zonoBundle object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(zB,type);
else
    [empty,res,S] = representsa_emptyObject(zB,type);
end
if empty; return; end

% dimension
n = dim(zB);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        res = representsa_(conZonotope(zB),'origin',tol);
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));
        
    case 'conPolyZono'
        res = true;
        if nargout == 2
            S = conPolyZono(zB);
        end

    case 'conZonotope'
        res = true;
        if nargout == 2
            S = conZonotope(zB);
        end

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));

    case 'halfspace'
        % polynomial zonotopes are bounded
        res = false;

    case 'interval'
        res = n == 1 || aux_isInterval(zB,tol);
        if nargout == 2 && res
            S = interval(zB);
        end

    case 'levelSet'
        res = true;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from zonoBundle to ' type ' not supported.']));
        end

    case 'polytope'
        res = true;
        if nargout == 2
            S = polytope(zB);
        end

    case 'polyZonotope'
        res = true;
        if nargout == 2
            S = polyZonotope(zB);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % obviously true
        res = true;
        if nargout == 2
            S = zB;
        end

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));

    case 'hyperplane'
        % zonotope bundles are bounded
        res = false;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonoBundle to ' type ' not supported.']));

    case 'emptySet'
        res = aux_emptySet(zB,tol);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % zonotope bundles are bounded
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isInterval(zB,tol)

    % one-dimensional zonoBundles are always intervals
    if dim(zB) == 1
        res = true; return
    end
    
    % empty sets can also be represented by intervals
    if representsa_(zB,'emptySet',tol)
        res = true; return
    end
    
    % all other cases: check individual zonotopes
    % (note: there are cases where the intersection is still an interval)
    res = true;
    for i=1:zB.parallelSets
        if ~representsa_(zB.Z{1},'interval',tol)
            res = false; return
        end
    end

end

function res = aux_emptySet(zB,tol)

    % full-empty zonotope bundle
    if zB.parallelSets == 0
        res = true; return
    end
    
    % check if any of the single zonotopes is empty
    for i = 1:length(zB.Z)
       if representsa_(zB.Z{i},'emptySet',eps)
           res = true;
           return;
       end
    end
    
    % check if the intersection of the zonotopes is empty
    res = representsa_(conZonotope(zB),'emptySet',tol);

end

% ------------------------------ END OF CODE ------------------------------
