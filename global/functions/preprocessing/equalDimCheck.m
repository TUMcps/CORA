function equalDimCheck(obj,S)
% equalDimCheck - checks if two objects (second argument can be class
%   array or numeric matrix) have the same dimension
%
% Syntax:  
%    equalDimCheck(obj,S)
%
% Inputs:
%    obj - contSet object
%    S   - contSet array or numeric matrix
%
% Outputs:
%    -
%
% Example:
%    obj1 = ellipsoid.generateRandom('Dimension',2);
%    obj2 = zonotope.generateRandom('Dimension',2);
%    obj3 = zonotope.generateRandom('Dimension',3);
%    equalDimCheck(obj1,[obj2,obj3]); % throws error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor Gassmann
% Written:      05-July-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% only check if macro is enabled
if CHECKS_ENABLED

    % check if dimensions match
    if isnumeric(S)
        if dim(obj) ~= size(S,1) && ~isscalar(S)
            throw(CORAerror('CORA:dimensionMismatch',obj,S));
        end
    
    elseif isa(obj,'interval') || isa(obj,'taylm') || isa(S,'interval') || isa(S,'taylm')
    
        if dim(obj) ~= dim(S)
            throw(CORAerror('CORA:dimensionMismatch',obj,S));
        end
    
    else
        ind_dim = dim(obj)==dim(S);
        if ~all(ind_dim)
            S_i = S(find(~ind_dim,1));
            throw(CORAerror('CORA:dimensionMismatch',obj,S_i));
        end
    
    end

end

%------------- END OF CODE --------------
