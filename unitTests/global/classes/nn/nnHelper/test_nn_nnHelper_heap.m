function res = test_nn_nnHelper_heap()
% test_nn_nnHelper_heap - tests nnHelper.heap
%
% Syntax:
%    res = test_nn_nnHelper_heap()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/heap

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

heapEntries = {5, 3, 2, 4, 7};
heap = nnHelper.heap(heapEntries);

v = 0;

while ~isempty(heap)
    vtemp = heap.pop().key;
    res = res && v < vtemp;
    v = vtemp;
end

end

% ------------------------------ END OF CODE ------------------------------
