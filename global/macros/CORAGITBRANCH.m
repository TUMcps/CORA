function gitbranch = CORAGITBRANCH()
% CORAGITBRANCH - returns the git branch of CORA
%
% Syntax:
%    gitbranch = CORAGITBRANCH()
%
% Inputs:
%    -
%
% Outputs:
%    gitbranch (string) - git branch of CORA
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Tobias Ladner
% Written:       22-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init dummy value
gitbranch = '';

% save current directory
currDir = pwd;

try
    % Proper Matlab-git interaction only available after R2023b
    if isMATLABReleaseOlderThan('R2024a')
        throw(CORAerror('CORA:specialError','This function is only available for Matlab version >=R2024a.'))
    end

    % switch to CORA root directory
    cd(CORAROOT);

    % make CORAROOT a safe directory to have access to the .git folder 
    % even though it does not belong to the user executing Matlab.
    matlab.git.addSafeDirectory(CORAROOT)
    
    % get git repo object
    repo = gitrepo;

    % read branch
    if isempty(repo.CurrentBranch)
       throw(MException('shared_cmlink:git:RepositoryNotFound','No git branch is available.')) % reusing identifier here.
    end
    gitbranch = repo.CurrentBranch.Name;

    % switch back to current directory
    cd(currDir)

catch ME
    % switch back to current directory
    cd(currDir)

    % check if error due to missing git repository
    if strcmp(ME.identifier,'shared_cmlink:git:RepositoryNotFound')
        % show warning
        CORAwarning('CORA:global','CORA does not have a git repository associated with its root folder.')
    else 
        % rethrow exception
        rethrow(ME)
    end
end

end

% ------------------------------ END OF CODE ------------------------------
