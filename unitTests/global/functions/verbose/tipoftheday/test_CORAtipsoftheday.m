function res = test_CORAtipsoftheday
% test_CORAtipsoftheday - unit test for test_CORAtipsoftheday
%
% Syntax:
%    res = test_CORAtipsoftheday
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       19-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get tips
tips = getCORAtipsoftheday();

% check if all tips can be displayed
for i=1:numel(tips)
    % build tip
    tip = tips{i};
    tip = compose(tip);
    tip = tip{1};

    % show tip
    fprintf('Tip #%i:\n',i)
    disp(tip)
    disp(' ')
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
