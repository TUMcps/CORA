function empty = writeSparseMatrixOptimized(fid,ind,var,tayMod)
% writeSparseMatrixOptimized - ???
%
% Syntax:
%    empty = writeSparseMatrixOptimized(fid,ind,var,tayMod)
%
% Inputs:
%    fid - identifier of the file to which the matrix is written
%    ind - ???
%    var - name of the matrix that is written
%    tayMod - ???
%
% Outputs:
%    empty - true/false whether index is empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: writeSparseMatrix, writeMatrix

% Authors:       ???
% Written:       ---
% Last update:   09-October-2024 (MW, add output argument, sprintf)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% variable
if contains(var,'Hf')
    out = 'outDyn';
elseif contains(var,'Hg')
    out = 'outCon';
else % for third-order tensors
    out = 'out';
end

empty = isempty(ind.row);
% loop over all non-empty entries
for i=1:length(ind.row)
    if tayMod
        str = sprintf('%s(%i,%i) = interval(%s(%i));',...
            var,ind.row(i),ind.col(i),out,ind.index(i));
    else
        str = sprintf('%s(%i,%i) = %s(%i);',...
            var,ind.row(i),ind.col(i),out,ind.index(i));
    end
    fprintf(fid, '%s\n', str);
end

% ------------------------------ END OF CODE ------------------------------
