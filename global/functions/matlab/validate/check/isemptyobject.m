function res = isemptyobject(S)
% isemptyobject - this function overloads the contSet/isemptyobject for
%    'numeric' variables
%
% Syntax:
%    res = isemptyobject(S)
%
% Inputs:
%    S - numeric vector/matrix
%
% Outputs:
%    res - true/false
%
% Example: 
%    x = [];
%    isemptyobject(x)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty, contSet/isemptyobject

% Authors:       Mark Wetzlinger
% Written:       03-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% numeric (point instead of set)
if isnumeric(S)
    if isempty(S)
        res = true;
    else
        res = false;
    end
else
    throw(CORAerror('CORA:notSupported',...
        'Function ''isemptyobject'' only supports contSet objects and numeric data types.'));
end

% ------------------------------ END OF CODE ------------------------------
