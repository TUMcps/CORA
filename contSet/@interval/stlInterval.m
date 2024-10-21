function res = stlInterval(I)
% stlInterval - converts an interval object to an stlInterval object
%
% Syntax:
%    res = stlInterval(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - stlInterval object
%
% Example:
%    I = interval(1,2);
%    res = stlInterval(1,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stlInterval

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% interval must be 1D
if ~all(size(I) == [1,1])
    throw(CORAerror('CORA:wrongValue',...
        'Interval must be one-dimensional.'));
end

% regular intervals are always closed
res = stlInterval(I.inf,I.sup,true,true);

% ------------------------------ END OF CODE ------------------------------
