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
res = all(size(arr) == [1,2]) && all(arr == [1,K]);

% maximum number of parts
arr = splitIntoNParts(K,K);
res(end+1,1) = all(size(arr) == [K,2]) && all(all(arr == [1:K; 1:K]'));

% even partition
N = 5;
arr = splitIntoNParts(K,N);
res(end+1,1) = all(size(arr) == [N,2]) ...
    && all(all(arr == [1 4; 5 8; 9 12; 13 16; 17 20]));

% uneven partition
K = 31;
N = 7;
arr = splitIntoNParts(K,N);
res(end+1,1) = all(size(arr) == [N,2]) && arr(1) == 1 && arr(end) == K;
res(end+1,1) = true;
for i=2:N
    % monotonically increasing
    if arr(i,2) < arr(i,1)
        res(end,1) = false; break
    end
    % no gaps
    if i > 1 && arr(i,1) ~= arr(i-1,2)+1
        res(end,1) = false; break
    end
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
