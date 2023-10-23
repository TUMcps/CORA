function equalDimCheck(S1,S2)
% equalDimCheck - checks if two objects are compatible, i.e., whether they
%    can execute a set-based operation
%
% Syntax:
%    equalDimCheck(S1,S2)
%
% Inputs:
%    S1 - contSet array or numeric matrix
%    S2 - contSet array or numeric matrix
%
% Outputs:
%    -
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
% See also: ---

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       05-July-2022
% Last update:   05-December-2022 (MW, fix matrix-case)
%                03-April-2023 (MW, integrate matrix-numeric case)
%                21-October-2023 (TL, cellfun)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only check if macro is enabled
if CHECKS_ENABLED

    if ismatrixset(S1) && isnumeric(S2)
        % operation between matrixset and matrix
        if ~all(dim(S1) == size(S2))
            throw(CORAerror('CORA:dimensionMismatch',S1,S2));
        end

    elseif (isnumeric(S1) || ismatrixset(S1)) && isa(S2,'contSet')
        % operation between matrix/matrix set and contSet
        % column dimension of matrix has to fit set dimension
        if size(S1,2) ~= dim(S2) && ~isscalar(S1)
            throw(CORAerror('CORA:dimensionMismatch',S1,S2));
        end

    elseif isnumeric(S2)
        % operation between set and matrix/vector/scalar
        % row dimension of matrix has to fit set dimension
        if dim(S1) ~= size(S2,1) && ~isscalar(S2)
            throw(CORAerror('CORA:dimensionMismatch',S1,S2));
        end
    
    elseif isa(S1,'interval') || isa(S1,'taylm') || isa(S1,'affine') || isa(S1,'zoo') ...
        || isa(S2,'interval') || isa(S2,'taylm') || isa(S2,'affine') || isa(S2,'zoo')
    
        % no class arrays for these classes
        if dim(S1) ~= dim(S2) && dim(S1) ~= 0 && dim(S2) ~= 0
            throw(CORAerror('CORA:dimensionMismatch',S1,S2));
        end

    elseif iscell(S1)
        cellfun(@(S1) equalDimCheck(S1,S2),S1);
    elseif iscell(S2)
        cellfun(@(S2) equalDimCheck(S1,S2),S2);
    
    else
        % operation between set and set
        ind_dim = dim(S1)==dim(S2);
        if ~all(ind_dim)
            S_i = S2(find(~ind_dim,1));
            throw(CORAerror('CORA:dimensionMismatch',S1,S_i));
        end
    
    end

end

% ------------------------------ END OF CODE ------------------------------
