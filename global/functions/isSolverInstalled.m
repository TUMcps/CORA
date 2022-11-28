function res = isSolverInstalled(name)
% isSolverInstalled - checks if specified solver is installed
%
% Syntax:  
%    res = isSolverInstalled(name)
%
% Inputs:
%    name - name of solver, currently available:
%           'mosek', 'sedumi', 'sdpt3', 'coneprog'
%
% Outputs:
%    res - true if installed, false otherwise
%
% Example: 
%    isSolverInstalled('mosek')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      12-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

switch lower(name)
    case 'mosek'
        res = logical(exist('mosekopt','file'));
    case 'sedumi'
        res = logical(exist('sedumi','file'));
    case 'sdpt3'
        res = logical(exist('sdpt3','file'));
    case 'coneprog'
        % is built-in starting R2020b
        res = exist('coneprog','file') > 0;
    otherwise
        res = false;
end

%------------- END OF CODE --------------