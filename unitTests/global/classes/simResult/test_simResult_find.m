function res = test_simResult_find
% test_simResult_find - unit test function for find
%
% Syntax:
%    res = test_simResult_find()
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
% Written:       22-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% init simResult with location
t = {[0; 1; 2],[2; 3; 4],[4; 5; 6],[6; 7; 8]};
x = {[0 0; 1 -1; 2 -2],[3 -2; 2 -1; 1 0],[1 0; 2 -1; 3 -1],[3 -1; 2 0; 1 1]};
loc = [1; 2; 1; 2];
simRes = simResult(x,t,loc);

% find all points in location 1
simRes_ = find(simRes,'location',1);

% check lengths
res(end+1,1) = length(simRes_) == 1 && length(simRes_.t) == 2 ...
    && length(simRes_.x) == 2 && length(simRes_.loc) == 2;

% check equality of time vector
res(end+1,1) = compareMatrices(simRes_.t{1},t{1}) ...
    && compareMatrices(simRes_.t{2},t{3});

% check equality of state vector
res(end+1,1) = compareMatrices(simRes_.x{1}',x{1}') ...
    && compareMatrices(simRes_.x{2}',x{3}');

% check equality of location vector
res(end+1,1) = all(cellfun(@(x) x == 1,simRes_.loc,'UniformOutput',true));

% same with location vector
t = {[0; 1; 2],[2; 3; 4],[4; 5; 6],[6; 7; 8]};
x = {[0 0; 1 -1; 2 -2],[3 -2; 2 -1; 1 0],[1 0; 2 -1; 3 -1],[3 -1; 2 0; 1 1]};
loc = [1 1; 1 2; 1 1; 1 2];
simRes = simResult(x,t,loc);

% find all points in location [1;2]
simRes_ = find(simRes,'location',[1;2]);

% check lengths
res(end+1,1) = length(simRes_) == 1 && length(simRes_.t) == 2 ...
    && length(simRes_.x) == 2 && length(simRes_.loc) == 2;

% check equality of time vector
res(end+1,1) = compareMatrices(simRes_.t{1},t{2}) ...
    && compareMatrices(simRes_.t{2},t{4});

% check equality of state vector
res(end+1,1) = compareMatrices(simRes_.x{1}',x{2}') ...
    && compareMatrices(simRes_.x{2}',x{4}');

% check equality of location vector
res(end+1,1) = all(cellfun(@(x) all(x' == [1;2]),simRes_.loc,'UniformOutput',true));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
