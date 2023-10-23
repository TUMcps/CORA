function res = test_simResult_extractHits
% test_simResult_extractHits - unit test function for extractHits
%
% Syntax:
%    res = test_simResult_extractHits()
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
% Written:       21-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty simResult
simRes = simResult();
[tHit,xHit,xHit_] = extractHits(simRes);
res = isempty(tHit) && isempty(xHit) && isempty(xHit_);

% simResult with trajectory (no specific location given)
t = {[0; 0.02; 0.05]};
x = {[1 1; 0.9 1.1; 0.8 1.2]};
simRes = simResult(x,t);
[tHit,xHit,xHit_] = extractHits(simRes);
res(end+1,1) = isempty(tHit) && isempty(xHit) && isempty(xHit_);

% simResult with one trajectory and given location
t = {[0; 0.02; 0.05]; [0.05; 0.08; 0.1]};
x = {[1 1; 0.9 1.1; 0.8 1.2]; [0.8 2.2; 0.7 2.3; 0.6 2.4]};
loc = [1; 2];
simRes = simResult(x,t,loc);

% no specific location
[tHit,xHit,xHit_] = extractHits(simRes);
res(end+1,1) = tHit == 0.05 ...
    && all(xHit{1} == [0.8; 1.2]) ...
    && all(xHit_{1} == [0.8; 2.2]);
% extract specific location before jump
[tHit,xHit,xHit_] = extractHits(simRes,2);
res(end+1,1) = isempty(tHit) && isempty(xHit) && isempty(xHit_);
% extract specific location after jump
[tHit,xHit,xHit_] = extractHits(simRes,[],2);
res(end+1,1) = tHit == 0.05 ...
    && all(xHit{1} == [0.8; 1.2]) ...
    && all(xHit_{1} == [0.8; 2.2]);

% simResult with location vector
t = {[0; 0.02; 0.05]; [0.05; 0.08; 0.1]};
x = {[1 1; 0.9 1.1; 0.8 1.2]; [0.8 2.2; 0.7 2.3; 0.6 2.4]};
loc = [1 1; 1 2];
simRes = simResult(x,t,loc);

% extract hits
[tHit,xHit,xHit_] = extractHits(simRes,[1;1]);
res(end+1,1) = tHit == 0.05 ...
    && all(xHit{1} == [0.8; 1.2]) ...
    && all(xHit_{1} == [0.8; 2.2]);

% wrong size of location vector
[tHit,xHit,xHit_] = extractHits(simRes,[1;1;1]);
res(end+1,1) = isempty(tHit) && isempty(xHit) && isempty(xHit_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
