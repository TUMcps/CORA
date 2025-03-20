function res = equalDimCheck(S1,S2,varargin)
% equalDimCheck - checks if two objects are compatible, i.e., whether they
%    can execute a set-based operation
%
% Syntax:
%    equalDimCheck(S1,S2)
%    equalDimCheck(S1,S2,returnValue)
%    res = equalDimCheck(S1,S2,returnValue)
%
% Inputs:
%    S1 - contSet array or numeric matrix
%    S2 - contSet array or numeric matrix
%    returnValue - true (return a value), false (error should be thrown)
%
% Outputs:
%    res - true/false
%
% Example:
%    M = [2 3; -1 -2];
%    C = capsule([1;1],[1;-1],0.5);
%    Z = zonotope([0;0;1],[1 -0.5 0; 0.4 0.6 1; 0 0.5 -1]);
% 
%    equalDimCheck(M,C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       05-July-2022
% Last update:   05-December-2022 (MW, fix matrix-case)
%                03-April-2023 (MW, integrate matrix-numeric case)
%                21-October-2023 (TL, cellfun)
%                30-September-2024 (MW, allow return value)
%                27-February-2025 (TL, allow broadcasting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

returnValue = setDefaultValues({false},varargin);

% one of the input arguments is a cell array
if iscell(S1)
    res = cellfun(@(S1) equalDimCheck(S1,S2,returnValue),S1,'UniformOutput',true);
    return
elseif iscell(S2)
    res = cellfun(@(S2) equalDimCheck(S1,S2,returnValue),S2,'UniformOutput',true);
    return
end

% only check if macro is enabled
res = aux_equalDimCheck(S1,S2,returnValue);
if ~returnValue && ~res
    throw(CORAerror('CORA:dimensionMismatch',S1,S2));
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_equalDimCheck(S1,S2,returnValue)

res = true;

if returnValue || CHECKS_ENABLED

    if isa(S1,'matrixSet') && isnumeric(S2)
        % operation between matrixset and matrix
        if ~all(dim(S1) == size(S2))
            res = false;
            return
        end

    elseif (isnumeric(S1) || isa(S1,'matrixSet')) && isa(S2,'contSet')
        % operation between matrix/matrix set and contSet
        % column dimension of matrix has to fit set dimension
        if size(S1,2) ~= dim(S2) && ~isscalar(S1)
            res = false;
            return
        end

    elseif isnumeric(S2)
        % operation between set and matrix/vector/scalar
        % row dimension of matrix has to fit set dimension
        if isscalar(dim(S1)) 
            % S1 is regular set
            if dim(S1) ~= size(S2,1) && ~isscalar(S2)
                res = false;
                return
            end
        else
            % S1 is matrix set, allow broadcasting
            dimsS1 = dim(S1);
            dimsS2 = size(S2,1:numel(dimsS1));
            if ~all(dimsS1 == dimsS2 | dimsS1 == 1 | dimsS2 == 1)
                res = false;
                return
            end            
        end
    
    elseif isa(S1,'interval') || isa(S1,'taylm') || isa(S1,'affine') || isa(S1,'zoo') ...
        || isa(S2,'interval') || isa(S2,'taylm') || isa(S2,'affine') || isa(S2,'zoo')
    
        % no class arrays for these classes
        if any(dim(S1) ~= dim(S2))
            res = false;
            return
        end
    
    else
        % operation between set and set
        ind_dim = dim(S1)==dim(S2);
        if ~all(ind_dim)
            % S_i = S2(find(~ind_dim,1));
            % throw(CORAerror('CORA:dimensionMismatch',S1,S_i));
            res = false;
            return
        end
    
    end

end

end

% ------------------------------ END OF CODE ------------------------------
