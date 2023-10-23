function res = test_batchCombinator()
% test_batchCombinator - unit test function for batchCombinator
%
% Syntax:
%    res = test_batchCombinator()
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
% See also: batchCombinator

% Authors:       Michael Eichelbeck
% Written:       02-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Test 1: 0 combinations -------------------------------------------

comb_state = struct;
[batch, comb_state] = batchCombinator(10, int16(13), 3, comb_state);
res = res & (isempty(batch) & comb_state.done==true);

% Test 2: 1 combination -------------------------------------------

comb_state = struct;
[batch, comb_state] = batchCombinator(10, int16(10), 10, comb_state);
expected = [1 2 3 4 5 6 7 8 9 10];
res = res & (all(batch==expected) & comb_state.done==true);

% Test 3: 10 combinations -------------------------------------------

comb_state = struct;

n = 5;
k = 3;

combs = zeros(10,3);
i = 1;

while true
    
    [batch, comb_state] = batchCombinator(n, int16(k), 3, comb_state);
    
    for j=1:length(batch(:,1))
        combs(i,:) = batch(j,:);
        i = i + 1;
    end

    if comb_state.done == true
        break;
    end

end

% Test if the correct number of combinations has been generated
% while their order is irrelvant

%are min and max indices correct
res_minmax = (max(max(combs)) == n);
res_minmax = res_minmax & (min(min(combs)) == 1);

%are all rows unique?
unique_rows = unique(combs, "rows");
res_unique = (size(unique_rows,1) == 10);

%does each row have only combinations without repetitions?
for i=1:size(combs, 1)
    unique_entries = unique(combs(i,:));
    res_unique = res_unique & (size(unique_entries,2) == k);
end

res = res & res_minmax & res_unique;

% ------------------------------ END OF CODE ------------------------------
