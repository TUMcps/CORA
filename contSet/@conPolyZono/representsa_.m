function [res,S] = representsa_(cPZ,type,tol,method,iter,splits)
% representsa_ - checks if a constrained polynomial zonotope can also be
%    represented by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(cPZ,type,tol,method,iter,splits)
%    [res,S] = representsa_(cPZ,type,tol,method,iter,splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', or 'all')
%    iter - number of iteration (integer > 0 or 'fixpoint')
%    splits - number of recursive splits (integer > 0)
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Example:
%    c = [0;0];
%    G = [1 0 1;0 1 1];
%    E = [1 0 2;0 1 1];
%    A = [1 -1 0; 0 -1 1];
%    b1 = [0; 1]; b2 = [0; 0];
%    EC = [2 0 1; 0 1 0];
%    cPZ1 = conPolyZono(c,G,E,A,b1,EC);
%    cPZ2 = conPolyZono(c,G,E,A,b2,EC);
%
%    res1 = representsa_(cPZ1,'emptySet',eps,'linearize',3,7)
%    res2 = representsa_(cPZ2,'emptySet',eps,'linearize',3,7)
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
    [empty,res] = representsa_emptyObject(cPZ,type);
else
    [empty,res,S] = representsa_emptyObject(cPZ,type);
end
if empty; return; end

% dimension
n = dim(cPZ);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));
        
    case 'conPolyZono'
        % obviously true
        res = true;
        if nargout == 2
            S = cPZ;
        end

    case 'conZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'halfspace'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'interval'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'polytope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'hyperplane'
        % constrained polynomial zonotopes cannot be unbounded (unless 1D,
        % where hyperplane is also bounded)
        res = n == 1;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conPolyZono to ' type ' not supported.']));

    case 'emptySet'
        res = aux_isEmptySet(cPZ,tol,method,splits,iter);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % constrained polynomial zonotopes cannot be unbounded
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isEmptySet(cPZ,tol,method,splits,iter)

    % check if independent generators are empty 
    if ~isempty(cPZ.GI)
       res = false; return; 
    end
    
    % check if constraints exist
    if isempty(cPZ.A)
       res = false; return; 
    end
        
    % try to contract the domain to the empty set -> set is empty
    temp = ones(length(cPZ.id),1);
    dom = interval(-temp,temp);
    
    D = contractPoly(-cPZ.b,cPZ.A,[],cPZ.EC,dom,method,iter,splits);
    
    res = representsa_(D,'emptySet',tol);

end

% ------------------------------ END OF CODE ------------------------------
