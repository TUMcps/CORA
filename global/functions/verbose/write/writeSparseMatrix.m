function empty = writeSparseMatrix(fid,M,var,varargin)
% writeSparseMatrix - write a sparse matrix to file
%
% Syntax:
%    empty = writeSparseMatrix(fid,M,var)
%    empty = writeSparseMatrix(fid,M,var,taylMod)
%
% Inputs:
%    fid - identifier of the file to which the matrix is written
%    M - symbolic matrix
%    var - name of the matrix that is written
%    taylMod - true/false whether inputs are Taylor models
%
% Outputs:
%    empty - true if matrix is empty, false otherwise
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: writeMatrix

% Authors:       Niklas Kochdumper
% Written:       15-July-2017
% Last update:   20-July-2017
%                24-January-2018
% Last revision: 09-October-2024 (MW, unify with writeSparseMatrixTaylorModel)

% ------------------------------ BEGIN CODE -------------------------------

% taylorModel or not?
taylMod = setDefaultValues({false},varargin);

[row,col] = find(M~=0);
empty = isempty(row);

% loop over all non-zero entries and print them one-by-one
if taylMod
    for i=1:length(row)
        fprintf(fid, '%s(%i,%i) = interval(%s);\n', ...
            var,row(i),col(i),bracketSubs(char(vpa(M(row(i),col(i))))));
    end
else
    for i=1:length(row)
        fprintf(fid,'%s(%i,%i) = %s;\n',...
            var,row(i),col(i),bracketSubs(char(M(row(i),col(i)))));
    end
end

% ------------------------------ END OF CODE ------------------------------
