function [empty,res,S_conv] = representsa_emptyObject(S,type)
% representsa_emptyObject - checks if a contSet class is a fully empty
%    object and returns an empty instance of class 'type'
%
% Syntax:
%    [empty,res] = representsa_emptyObject(S,type)
%    [empty,res,S_conv] = representsa_emptyObject(S,type)
%
% Inputs:
%    S - contSet object
%    type - contSet class
%
% Outputs:
%    empty - true/false whether S is fully empty
%    res - true/false whether S can be represented by 'type'
%    S_conv - empty object of class type
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if set is a fully empty object
empty = isemptyobject(S);

switch type
    case {'origin','point','fullspace','parallelotope'}
        res = false;
        S_conv = [];
        
    case 'hyperplane'
        res = true;
        if nargout == 3
            S_conv = contSet.initEmptySet('conHyperplane');
        end

    otherwise
        res = true;
        if nargout == 3
            S_conv = contSet.initEmptySet(type);
        end
end 

% ------------------------------ END OF CODE ------------------------------
