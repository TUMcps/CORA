function res = isequal(I,S,varargin)
% isequal - checks if an interval is equal to another set or point
%
% Syntax:
%    res = isequal(I,S)
%    res = isequal(I,S,tol)
%
% Inputs:
%    I - interval object
%    S - contSet object, numeric
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

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({eps},varargin);

% input argument check
inputArgsCheck({{I,'att',{'interval','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% ensure that numeric is second input argument
[I,S] = reorderNumeric(I,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < I.precedence
    res = isequal(S,I,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(I,S,true)
    res = false;
    return
end

% interval-interval case
if isa(S,'interval')
    res = aux_isequal_interval(I,S,tol);
    return
end

% interval-vector case
if isnumeric(S)
    S = interval(S);
    res = aux_isequal_interval(I,S,tol);
    return
end

throw(CORAerror('CORA:noops',I,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_interval(I,S,tol)

% assume false
res = false;

% check for equal dimensions
if any(dim(I) ~= dim(S))
    return
end

% read infima and suprema
lb1 = infimum(I); ub1 = supremum(I);
lb2 = infimum(S); ub2 = supremum(S);

% indices with infinity values
idxInf_lb1 = isinf(lb1); idxInf_ub1 = isinf(ub1);
idxInf_lb2 = isinf(lb2); idxInf_ub2 = isinf(ub2);

% checks
if ~all(idxInf_lb1 == idxInf_lb2,"all") || ~all(idxInf_ub1 == idxInf_ub2,"all")
    % if same entries are minus/plus infinity
    return
elseif any(lb1(idxInf_lb1) ~= lb2(idxInf_lb2),"all") || any(ub1(idxInf_ub1) ~= ub2(idxInf_ub2),"all")
    % if infinity entries are equal
    return
elseif ~all(withinTol(lb1(~idxInf_lb1),lb2(~idxInf_lb2),tol),"all") ...
        || ~all(withinTol(ub1(~idxInf_ub1),ub2(~idxInf_ub2),tol),"all")
    % if other values are equal
    return
end

% assume true
res = true;

end

% ------------------------------ END OF CODE ------------------------------
