function [Enew,Gnew] = removeRedundantExponents(E,G)
% removeRedundantExponents - add up all generators that belong to terms
%    with identical exponents
%
% Syntax:
%    [Enew, Gnew] = removeRedundantExponents(E,G)
%
% Inputs:
%    E - matrix containing the exponent vectors
%    G - generator matrix
%
% Outputs:
%    Enew - modified exponent matrix
%    Gnew - modified generator matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       25-June-2018
% Last update:   21-April-2020 (remove zero-length generators)
%                23-June-2022 (performance optimizations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% return directly if G is empty
if isempty(G)
    Enew = E;
    Gnew = G;
    return;
end

% remove zero-length generators
idxD = any(G,1);

% skip if all non-zero
if ~all(idxD)
    if all(~idxD)
        Enew = zeros(size(E,1),1);
        Gnew = zeros(size(G,1),1);
        return;
    else
        G = G(:,idxD);
        E = E(:,idxD);
    end
end

% create a deterministic random hash vector
rs = RandStream('mt19937ar', 'Seed', 0); % to not interfer with the outside
hashVec = rs.rand(1, size(E, 1));
hashMat = (hashVec*E)';

% sort the exponent vectors according to the hash value
[hashes, ind] = sortrows(hashMat);
ind = ind'; % row

% test if all (sorted) hashes are different
uniqueHashIdx = [hashes(1:end-1) ~= hashes(2:end); true] ... % succ 
    & [true; hashes(1:end-1) ~= hashes(2:end)]; % pred

% if so, return directly
numUnique = sum(uniqueHashIdx);
if numUnique == size(E, 2)
    Enew = E;
    Gnew = G;
    return
end

% initialize new matrices
Enew = zeros(size(E));
Gnew = zeros(size(G));

% copy unique hashes
uniqueColumns = ind(uniqueHashIdx);
Enew(:, 1:numUnique) = E(:, uniqueColumns);
Gnew(:, 1:numUnique) = G(:, uniqueColumns);
current = numUnique + 1;

% continue with potential redundancies
ind = ind(~uniqueHashIdx);

E = E(:, ind);
G = G(:, ind);
hashMat = hashMat(ind);
ind = 1:length(ind);

% first entry
ind_i = ind(1);
exp_c = E(:,ind_i);
hash_c = hashMat(ind_i);
Enew(:,current) = exp_c;
Gnew(:,current) = G(:,ind_i);

% loop over all exponent vectors
for ind_i = ind(2:end)
    hash_i = hashMat(ind_i);
    exp_i = E(:,ind_i);

    if hash_c == hash_i && all(exp_c == exp_i)           
        Gnew(:,current) = Gnew(:,current) + G(:,ind_i);
    else
        current = current + 1;
        exp_c = exp_i;
        hash_c = hash_i;
        Enew(:,current) = exp_c;
        Gnew(:,current) = G(:,ind_i);
    end
end 


% truncate exponent and generator matrix
Enew(:,current+1:end) = [];
Gnew(:,current+1:end) = [];

% ------------------------------ END OF CODE ------------------------------
