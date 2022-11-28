function res = isequal(I1,I2,varargin)
% isequal - checks if two intervals are equal
%
% Syntax:  
%    res = isequal(I1,I2,tol)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([1; -1; 0], [4; 2; 1]);
%    I2 = interval([1; 0; 0], [3.5; 2; 1]);
%    res = isequal(I1,I2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  12-March-2021 (MW, add dimension mismatch)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
tol = setDefaultValues({eps},varargin{:});

% input argument check
inputArgsCheck({{I1,'att','interval'};
                {I2,'att','interval'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

if dim(I1) ~= dim(I2)
    throw(CORAerror('CORA:dimensionMismatch',I1,I2));    
else
    res = all(abs(infimum(I1) - infimum(I2)) < tol) && ... % infima
        all(abs(supremum(I1) - supremum(I2)) < tol); % suprema
end

%------------- END OF CODE --------------