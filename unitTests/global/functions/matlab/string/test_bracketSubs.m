function res = test_bracketSubs
% test_bracketSubs - unit test function for substitution of brackets with
%    letters 'L' and 'R'
%
% Syntax:
%    res = test_bracketSubs()
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

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init state string
str = 'xL1R';
str_ = bracketSubs(str);
str_true = 'x(1)';
res = strcmp(str_,str_true);

% input string
str = 'uL16R';
str_ = bracketSubs(str);
str_true = 'u(16)';
res(end+1,1) = strcmp(str_,str_true);

% longer name
str = 'abcdeL99R';
str_ = bracketSubs(str);
str_true = 'abcde(99)';
res(end+1,1) = strcmp(str_,str_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
