function res = releaseCORA()
% releaseCORA - create a new CORA release
%
% Syntax:
%    releaseCORA()
%
% Inputs:
%    -
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Tobias Ladner
% Written:       09-February-2024
% Last update:   07-March-2024 (TL, updates)
%                17-October-2025 (TL, rewrite devgeneral -> PUBLIC)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp(' ')
disp('-------------------------------------------------------------------')
disp('                           CORA RELEASE                            ')
disp('-------------------------------------------------------------------')
disp(' ')

% change dir
currentdir = pwd;
cd(CORAROOT);

% make sure all branches are up-to-date
disp('-------------------------------------------------------------------')
disp("Checking if all branches are up-to-date:")
branches = {'PUBLIC','devgeneral','master'};
aux_checkBranchesUpToDate(branches)

% read commit message
disp('-------------------------------------------------------------------')
aux_switchBranch('devgeneral')
changelogfile = sprintf('./global/functions/release/releaseNotes.txt');
newCORAVERSION = aux_readChangelog(changelogfile);

% make pre-release commit on devgeneral
disp('-------------------------------------------------------------------')
disp("Preparing pre-release on devgeneral:")
aux_switchBranch('devgeneral')
disp("Updating CORAVERSION:")
aux_updateCORAVERSION(newCORAVERSION);
system(sprintf('git add %s',which('CORAVERSION')));
system(sprintf('git add %s',changelogfile));
system(sprintf('git commit -m "pre-release for %s"',newCORAVERSION));
aux_awaitUser();

% merge devgeneral into PUBLIC
disp('-------------------------------------------------------------------')
disp("Merging devgeneral into PUBLIC:")
aux_switchBranch('PUBLIC')
system('git merge --squash devgeneral');

% delete manual on PUBLIC
aux_removeManual()

% commit release
disp('-------------------------------------------------------------------')
disp('Commit release on PUBLIC:')
system(sprintf('git commit -F %s',changelogfile));
aux_awaitUser();

% tag version
disp('-------------------------------------------------------------------')
disp("Tagging..")
tagname = newCORAVERSION(6:end);
system(sprintf('git tag %s', tagname));

% merge with all branches
disp('-------------------------------------------------------------------')
disp("Merging branches into each other:")
aux_awaitUser();
aux_mergeBranches('PUBLIC', 'devgeneral','-s ours')
aux_mergeBranches('devgeneral', 'master')
disp("- all branches merged!")

% push branches
disp('-------------------------------------------------------------------')
disp("Pushing branches to GitLab:")
aux_awaitUser();
aux_pushBranch('origin', 'PUBLIC')
aux_pushBranch('origin', 'devgeneral')
aux_pushBranch('origin', 'master')
system(sprintf('git push origin %s', tagname));

disp('-------------------------------------------------------------------')
disp("Pushing branches to GitHub:")
aux_awaitUser();
aux_pushBranch('github', 'PUBLIC', 'master')
system(sprintf('git push github %s', tagname));

disp('-------------------------------------------------------------------')
disp("Send message in slack:")
aux_slackmessage(changelogfile)
aux_awaitUser();

% reset directory
cd(currentdir);

disp('-------------------------------------------------------------------')
fprintf("'%s' release successful!\n", CORAVERSION)
disp('-------------------------------------------------------------------')

res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_checkBranchesUpToDate(branches)
    for i=1:numel(branches)
        branch = branches{i};
        fprintf("- Checking branch '%s'\n", branch);
        aux_switchBranch(branch)
        system('git pull');
        system('git lfs pull');
    end
end

function aux_switchBranch(branch)
    try
        res = system(sprintf('git checkout %s', branch));
        if res ~= 0
            throw(CORAerror('CORA:specialError',sprintf('Unable to checkout branch ''%s''. See above.', branch)))
        end
    catch ME
        if strcmp(ME.identifier, 'CORA:specialError') || ...
           ~contains(ME.cause{1}.message,'is already checked out')
            rethrow(ME);
        end
    end
end

function [newCORAVERSION] = aux_readChangelog(changelogfile)
    % prepare commit message
    disp('Prepare commit message:')
    disp('- Gathering commit history on devgeneral.. (edit subsequent command if not sufficient and run in CORAROOT)')
    commithistoryfile = strrep(changelogfile,'releaseNotes','commitHistory');
    commithistorycommand = sprintf('git log -100 --pretty=oneline --format=%%B > %s', commithistoryfile);
    fprintf('\n  %s\n\n', commithistorycommand);
    system(commithistorycommand);
    fprintf('- Commit history available in: %s\n', commithistoryfile)
    fprintf('- Write final release commit message into: %s\n', changelogfile)
    % open them in reverse to show commit history file to user
    fclose(fopen(changelogfile,'a')); % create file if not already
    edit(changelogfile,commithistoryfile)
    aux_awaitUser()

    % read commit message
    disp("Reading commit message..")
    commitmsg = fileread(changelogfile);
    lines = splitlines(commitmsg);

    % read new CORAVERSION
    newCORAVERSION = lines{1};

    % check new CORAVERSION
    fprintf("- Current CORAVERSION: '%s'\n", CORAVERSION)
    fprintf("- New CORAVERSION:     '%s'\n", newCORAVERSION)
    aux_awaitUser()

    % sort changelog
    if ~isempty(lines{end})
        lines{end+1} = '';
    end
    lines(3:end-1) = sort(lines(3:end-1));
    commitmsg = strjoin(lines,newline);

    % check changelog
    disp('Commit message:')
    disp('---')
    disp(commitmsg)
    disp('---')
    aux_awaitUser()

    % write back to file
    fileID = fopen(changelogfile,'w');
    fprintf(fileID,commitmsg);
    fclose(fileID);

end

function aux_awaitUser()
    answer = input('Continue? [y/n] >> ','s');
    if ~isequal(answer,'y')
        throw(CORAerror('CORA:specialError','Aborted.'))
    end
end

function aux_updateCORAVERSION(newCORAVERSION)
    % read file
    fileCORAVERSION = which("CORAVERSION");
    text = fileread(fileCORAVERSION);

    % replace CORAVERSION
    text = strrep(text, CORAVERSION,newCORAVERSION);

    % split lines 
    lines = splitlines(text);
    
    % update last update line
    today = datetime;
    today.Format = 'dd-MMMM-yyyy';
    lines{21} = ['% Last update:   ' char(today)];

    % write back to file
    fileID = fopen(fileCORAVERSION,'w');
    % for some reason Matlab doesn't want it to be written at once ...
    for i=1:numel(lines)-1
        fprintf(fileID, '%s\n',lines{i});
    end
    fclose(fileID);
end

function aux_removeManual()
    disp('Removing manual from PUBLIC:')
    manualfolder = sprintf('%s/manual', CORAROOT);
    rmpath(genpath(manualfolder))
    rmdir(manualfolder,'s');
    system('git add manual')
end

function aux_mergeBranches(branchFrom, branchTo, parameters)
    if nargin < 3
        parameters = '--no-ff';
    end
    fprintf('- Merging branch ''%s'' into ''%s'' (using ''%s'')\n', branchFrom, branchTo, parameters)
    aux_switchBranch(branchTo)
    system(sprintf('git merge %s %s', branchFrom, parameters));
    aux_awaitUser();
end

function aux_pushBranch(remote, branch, branch_remote)
    if nargin < 3
        branch_remote = branch;
    end

    fprintf('- Pushing ''%s'' to remote ''%s'' with branch name ''%s'':\n', branch, remote, branch_remote);
    system(sprintf('git push %s %s:%s', remote, branch, branch_remote));
end

function aux_slackmessage(changelogfile)

    % read commit message
    commitmsg = fileread(changelogfile);
    
    % add 'now released' text
    lines = splitlines(commitmsg);
    lines{1} = sprintf('%s %s', lines{1}, 'now released!');
    commitmsg = strjoin(lines,newline);

    % check changelog
    disp('---')
    disp(commitmsg)
    disp('Visit https://cora.in.tum.de/ for a surprise.')
    disp('---')
end

% ------------------------------ END OF CODE ------------------------------
