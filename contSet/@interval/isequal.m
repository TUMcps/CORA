function res = isequal(I1,I2,varargin)
% isequal - checks if two intervals are equal
%
% Syntax:
%    res = isequal(I1,I2)
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

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   12-March-2021 (MW, add dimension mismatch)
%                03-December-2022 (MW, add check for infinity)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% input argument check
inputArgsCheck({{I1,'att','interval'};
                {I2,'att','interval'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% assume false
res = false;

% check for equal dimensions
if any(dim(I1) ~= dim(I2))
    return
end

% read infima and suprema
lb1 = infimum(I1); ub1 = supremum(I1);
lb2 = infimum(I2); ub2 = supremum(I2);

% indices with infinity values
idxInf_lb1 = isinf(lb1); idxInf_ub1 = isinf(ub1);
idxInf_lb2 = isinf(lb2); idxInf_ub2 = isinf(ub2);

% checks
if ~all(all(idxInf_lb1 == idxInf_lb2)) || ~all(all(idxInf_ub1 == idxInf_ub2))
    % if same entries are minus/plus infinity
    return
elseif any(lb1(idxInf_lb1) ~= lb2(idxInf_lb2)) || any(ub1(idxInf_ub1) ~= ub2(idxInf_ub2))
    % if infinity entries are equal
    return
elseif ~all(withinTol(lb1(~idxInf_lb1),lb2(~idxInf_lb2),tol)) ...
        || ~all(withinTol(ub1(~idxInf_ub1),ub2(~idxInf_ub2),tol))
    % if other values are equal
    return
end

% assume true
res = true;

% ------------------------------ END OF CODE ------------------------------
