function res = isemptyobject(I)
% isemptyobject - checks whether an STL interval is empty
%
% Syntax:
%    res = isemptyobject(I)
%
% Inputs:
%    I - stlInterval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I = stlInterval(1,2);
%    isemptyobject(I); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isempty(I.lower);

% ------------------------------ END OF CODE ------------------------------
