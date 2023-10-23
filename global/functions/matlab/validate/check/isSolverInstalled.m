function res = isSolverInstalled(varargin)
% isSolverInstalled - checks if any of the specified solvers is installed
%
% Syntax:
%    res = isSolverInstalled(name1,name2,...)
%
% Inputs:
%    name(.) - name of solver, currently available:
%           'mosek', 'sedumi', 'sdpt3', 'coneprog', 'gurobi'
%
% Outputs:
%    res - true any of the specified is installed, false otherwise
%
% Example: 
%    isSolverInstalled('mosek','gurobi')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       12-May-2022
% Last update:   17-April-2023 (VG, added gurobi, check if any given solvers is supported)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i=1:length(varargin)
    name = varargin{i};
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
        case 'gurobi'
            res = logical(exist('gurobi','file'));
        otherwise
            res = false;
    end
    if res
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
