function res = compareIntervals(as,bs)
% compareIntervals - helper function for comparing stlInterval arrays
%
% Syntax:
%    res = compareIntervals(as,bs)
%
% Inputs:
%    as - array of stlInterval
%    bs - array of stlInterval
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       20-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isequal(size(as),size(bs))
    res = false;
    return
end

res = all(arrayfun(@(a,b) a == b,as,bs));

% ------------------------------ END OF CODE ------------------------------
