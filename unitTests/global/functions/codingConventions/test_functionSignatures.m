function res = test_functionSignatures()
% test_functionSignatures - tests if all functions signatures are valid
%
% Syntax:
%    res = test_functionSignatures()
%
% Inputs:
%    -
%
% Outputs:
%    res - whether all files are valid
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: functionSignatures.json

% Authors:       Tobias Ladner
% Written:       05-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% save current directory
currDir = pwd;

% search for functionSignagures.json
files = findfiles(CORAROOT,true,'','json');
files = files(arrayfun(@(file) strcmp(file.name,'functionSignatures.json'),files));

% validate all function signatures
try
    for i=1:numel(files)
        % check file i
        file = files(i);
        % functionSignatures.json has to be in the current folder ...
        cd(file.folder)
        % call validate function
        issues = validateFunctionSignaturesJSON;
        % assert no issue
        assertLoop(isempty(issues),[file.folder filesep file.name '.'],[],i);
    end
catch ME
    cd(currDir);
    rethrow(ME)
end

% switch back to current directory
cd(currDir);

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
