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
%    equalDimCheck(C,Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      05-July-2022
% Last update:  05-December-2022 (MW, fix matrix-case)
% Last revision:---

%------------- BEGIN CODE --------------

% only check if macro is enabled
if CHECKS_ENABLED

    % check if dimensions match
    if isnumeric(S1) || ismatrixset(S1)
        % operation between matrix/vector/scalar and set
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
        if dim(S1) ~= dim(S2)
            throw(CORAerror('CORA:dimensionMismatch',S1,S2));
        end
    
    else
        % operation between set and set
        ind_dim = dim(S1)==dim(S2);
        if ~all(ind_dim)
            S_i = S2(find(~ind_dim,1));
            throw(CORAerror('CORA:dimensionMismatch',S1,S_i));
        end
    
    end

end

%------------- END OF CODE --------------
