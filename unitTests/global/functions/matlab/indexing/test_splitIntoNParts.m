function res = test_splitIntoNParts
% test_splitIntoNParts - unit test function for splitIntoNParts
%
% Syntax:
%    res = test_splitIntoNParts()
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
% Written:       19-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% minimum number of parts
K = 20;
arr = splitIntoNParts(K,1);
assert(all(size(arr) == [1,2]) && all(arr == [1,K]));

% maximum number of parts
arr = splitIntoNParts(K,K);
assert(all(size(arr) == [K,2]) && all(all(arr == [1:K; 1:K]')));

% even partition
N = 5;
arr = splitIntoNParts(K,N);
assert(all(size(arr) == [N,2]) ...
    && all(arr == [1 4; 5 8; 9 12; 13 16; 17 20],"all"));

% uneven partition
K = 31;
N = 7;
arr = splitIntoNParts(K,N);
assert(all(size(arr) == [N,2]) && arr(1) == 1 && arr(end) == K);

for i=2:N
    % monotonically increasing
    assertLoop(arr(i,1) < arr(i,2),i)

    % no gaps
    assertLoop(i > 1 && arr(i-1,1)+1 ~= arr(i,2),i)
end

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
