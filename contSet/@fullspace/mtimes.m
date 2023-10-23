function res = mtimes(factor1,factor2)
% mtimes - overloaded '*' operator for the linear map of a full-dimensional
%    space
%    case R^0: can only be multiplied by 0 (not representable in MATLAB)
%
% Syntax:
%    fs = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - fullspace object, numerical scalar/matrix
%    factor2 - fullspace object, numerical scalar/matrix
%
% Outputs:
%    res - linearly mapped full-dimensional space
%
% Example: 
%    fs = fullspace(2);
%    M = [2 1; -1 3];
%    M*fs
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check dimensions
equalDimCheck(factor1,factor2);

% find the fullspace object
[fs,M] = findClassArg(factor1,factor2,'fullspace');

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported','Linear map of R^0 not supported'));
end

if isscalar(M)
    % multiplication with scalar
    if withinTol(M,0)
        res = interval(zeros(fs.dimension,1),zeros(fs.dimension,1));
    else
        % all other scalar values: keep fs as is
        res = fs;
    end

elseif isnumeric(M)
    zerorows = ~any(M,2);
    if any(zerorows)
        res = interval(-Inf(size(M,1),1),Inf(size(M,1),1));
        res(zerorows) = 0;
    elseif size(M,1) ~= fs.dimension
        % all other matrices: still fullspace
        res = fullspace(size(M,1));
    else
        % square matrix
        res = fs;
    end

elseif isa(M,'intervalMatrix') || isa(M,'matZonotope') || isa(M,'matPolytope')
    % throw error for now...
    throw(CORAerror('CORA:noops',M,fs));   

end

% ------------------------------ END OF CODE ------------------------------
