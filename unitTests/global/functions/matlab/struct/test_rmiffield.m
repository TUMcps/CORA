function res = test_rmiffield
% test_rmiffield - unit test function for rmiffield
%
% Syntax:
%    res = test_rmiffield()
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
% Written:       17-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty case
S = struct();
label = 'test';
S = rmiffield(S,label);
resvec(end+1) = ~isfield(S,label);

% normal case
S = struct();
label1 = 'label1';
label2 = 'label2';
S.(label1) = 'value1';
S.(label2) = 'value2';
resvec(end+1) = isfield(S,label1);
resvec(end+1) = isfield(S,label2);
S = rmiffield(S,label1);
resvec(end+1) = ~isfield(S,label1);
resvec(end+1) = isfield(S,label2);

% another call should not change the result
S = rmiffield(S,label1);
resvec(end+1) = ~isfield(S,label1);
resvec(end+1) = isfield(S,label2);

% combine results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
