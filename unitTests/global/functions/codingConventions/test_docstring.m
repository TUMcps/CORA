function res = test_docstring()
% test_docstring - legacy, moved to test_codingConventions
%
% Syntax:
%    res = test_docstring()
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
% See also: test_codingConventions

% Authors:       Tobias Ladner
% Written:       18-November-2022
% Last update:   05-February-2025 (moved to test_codingConventions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% show CORAwarning
CORAwarning('CORA:deprecated','function','test_docstring','CORA v2025.1.0','Please use test_codingConventions instead.','This change was made to desribe the purpose of this function more clearly. Please rerun the test using test_codingConventions.')

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
