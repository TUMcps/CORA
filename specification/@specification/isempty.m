function res = isempty(spec)
% isempty - checks if a specification object is empty
%
% Syntax:
%    res = isempty(spec)
%
% Inputs:
%    spec - specification object
%
% Outputs:
%    res - true/false
%
% Example:
%    spec = specification();
%    res = isempty(spec)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       02-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(spec) == 0
    res = true; return
end

res = any(arrayfun(@(x) isnumeric(x.set) && isempty(x.set),...
    spec,'UniformOutput',true));

% ------------------------------ END OF CODE ------------------------------
