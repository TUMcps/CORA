function [res,S] = representsa_(ls,type,tol,varargin)
% representsa_ - checks if a level set can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(ls,type,tol)
%    [res,S] = representsa_(ls,type,tol)
%
% Inputs:
%    ls - levelSet object
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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       19-July-2023
% Last update:   17-January-2023 (TL, emptySet & fullspace)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(ls,type);
else
    [empty,res,S] = representsa_emptyObject(ls,type);
end
if empty; return; end

% dimension
n = dim(ls);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));
        
    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'conZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'halfspace'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'interval'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'levelSet'
        % obviously true
        res = true;
        if nargout == 2
            S = ls;
        end

    case 'polytope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'hyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'emptySet'
        res = aux_representsa_emptySet(ls);

    case 'fullspace'
        res = aux_representsa_fullspace(ls);

end

end


% Auxiliary functions -----------------------------------------------------


function res = aux_representsa_emptySet(ls) 
    % assume true
    res = true;

    % evaluate equations
    n = dim(ls);
    I = interval.Inf(n);
    O = ls.funHan(I);
    O = interval(O); % in case its just numeric

    % convert compare operator to cell
    compOp = ls.compOp;
    if ~iscell(compOp)
        compOp = {compOp};
    end

    % check comparator
    for i=1:numel(compOp)
        switch compOp{i}
            case '=='
                if contains(O(i),0)
                    % possibly not empty
                    res = false;
                end

            case '<='
                if O.inf(i) <= 0
                    % possibly not empty
                    res = false;
                end

            case '<'
                if O.inf(i) < 0
                    % possibly not empty
                    res = false;
                end


            otherwise
                throw(CORAerror('CORA:specialError',sprintf("Unknown compare operator: '%s'.",compOp{i})));
        end
    end
end

function res = aux_representsa_fullspace(ls) 
    
    % assume true
    res = true;

    % evaluate equations
    n = dim(ls);
    I = interval.Inf(n);
    O = ls.funHan(I);
    O = interval(O); % in case its just numeric

    % convert compare operator to cell
    compOp = ls.compOp;
    if ~iscell(compOp)
        compOp = {compOp};
    end

    % check comparator
    for i=1:numel(compOp)
        switch compOp{i}
            case '=='
                if ~isequal(O(i),interval(0))
                    % possibly not fullspace
                    res = false;
                end

            case '<='
                if O.sup(i) > 0
                    % possibly not fullspace
                    res = false;
                end

            case '<'
                if O.sup(i) >= 0
                    % possibly not fullspace
                    res = false;
                end


            otherwise
                throw(CORAerror('CORA:specialError',sprintf("Unknown compare operator: '%s'.",compOp{i})));
        end
    end
end


% ------------------------------ END OF CODE ------------------------------
