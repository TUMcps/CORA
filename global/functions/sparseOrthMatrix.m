function Q = sparseOrthMatrix(d)
% sparseOrthMatrix - generates a sparse orthogonal matrix
%
% Syntax:  
%    Q = sparseOrthMatrix(d)
%
% Inputs:
%    d - dimension
%
% Outputs:
%    Q - sparse orthogonal matrix of dimension dim
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      07-Oct-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% generate blocks
min_block = 2;
max_block = 4;
remaining = d;
i = 1;
while true
    if i == 1
        % fix one block to size 2
        if d <= 4
            block(i) = d;
            break
        end
        next = 2;
    else
        next = min_block + floor(rand(1)*(max_block - min_block + 1));
    end
    if remaining - next == 0 || ...
            (remaining >= next && remaining - next >= min_block)
        block(i) = next;
        i = i + 1;
        remaining = remaining - next;
    end
    if remaining == 0
        break
    end
end

% write Q in block structure
blocks = length(block);
Q = zeros(d);
curr_start = 1;
curr_end = 0;
for i=1:blocks
    [Qblock,~] = qr(2*rand(block(i)) - 1); % took 2*rand - 1 from randomLinSys
    curr_end = curr_end + block(i);
    Q(curr_start:curr_end,curr_start:curr_end) = Qblock;
    curr_start = curr_start + block(i);
end

% reorder columns and rows -> no more strict block structure
reOrderCols = randperm(d);
reOrderRows = randperm(d);
Q = Q(reOrderRows,:);
Q = Q(:,reOrderCols);

% reorder such that last row has two non-zero entries
if nnz(Q(end,:)) ~= 2
    order = 1:length(Q);
    for i=1:length(Q)-1
        if nnz(Q(i,:)) == 2
            order(i) = length(Q);
            order(end) = i;
            Q = Q(order,:);
            break
        end
    end
end

% check if Q really orthogonal
if ~(abs(abs(det(Q)) - 1) < 1e-9 && all(abs(vecnorm(Q)-1) < 1e-9))
    error("Matrix not orthogonal!");
end



end


%------------- END OF CODE --------------