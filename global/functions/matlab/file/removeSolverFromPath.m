function [res,path2solver] = removeSolverFromPath(name)
% removeSolverFromPath - returns the CORA root path
%
% Syntax:
%    removeSolverFromPath
%
% Inputs:
%    name - name of solver, currently available:
%           'mosek'
%
% Outputs:
%    res - true/false whether solver could be removed
%    path2solver - cell-array of path to all folders which were removed
%
% Example:
%    [res,path2solver] = removeSolverFromPath('mosek')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isSolverInstalled.m

% Authors:       Mark Wetzlinger
% Written:       14-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init output arguments
res = true;
path2solver = {};

switch lower(name)
    case 'mosek'
        % path to solver
        temp = which('mosekdiag');
        while ~isempty(temp)
            % mosekdiag can be in ../toolbox/r2017a and ../toolbox/r2017aom
            try
                % directory
                temp = fileparts(temp);
                % remove from path
                rmpath(temp);
                res = true;
                % add to list of removed folders
                path2solver{end+1,1} = temp;
            catch ME
                % some issue...
                res = false; return
            end
            % check if another instance of mosekdiag is on the path
            temp = which('mosekdiag');
        end

    otherwise
        res = false;
end

% ------------------------------ END OF CODE ------------------------------
