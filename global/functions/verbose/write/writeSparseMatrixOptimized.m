function writeSparseMatrixOptimized(ind,var,fid,tayMod)
% writeSparseMatrixOptimized - ???
%
% Syntax:
%    writeSparseMatrixOptimized(ind,var,fid,tayMod)
%
% Inputs:
%    ind - ???
%    var - name of the matrix that is written
%    fid - identifier of the file to which the matrix is written
%    tayMod - ???
%
% Outputs:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       ???
% Written:       ---
% Last update:   ---
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

% loop over all non-empty entries
for i = 1:length(ind.row)
    if ~tayMod
        str = [var,'(',num2str(ind.row(i)),',',num2str(ind.col(i)), ...
               ') = ',out,'(',num2str(ind.index(i)),');'];
    else
        str = [var,'(',num2str(ind.row(i)),',',num2str(ind.col(i)), ...
               ') = interval(',out,'(',num2str(ind.index(i)),'));'];
    end
      
    fprintf(fid, '%s\n', str);
end

% ------------------------------ END OF CODE ------------------------------
