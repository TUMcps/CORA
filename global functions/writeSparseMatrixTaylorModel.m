function writeSparseMatrixTaylorModel(M,var,fid)
% writeSparseMatrixTaylorModel - write a sparse matrix in such a way that
%                                the inputs are taylor models and the
%                                output is a interval matrix
%
% Syntax:  
%    writeSparseMatrixTaylorModel(M,var,fid)
%
% Inputs:
%    M - symbolic matrix
%    var - name of the matrix that is written
%    fid - identifier of the file to which the matrix is written
%
% Outputs:
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      15-July-2017
% Last update:  20-July-2017
%               24-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

%write each row
[row,col] = find(M~=0);

for i=1:length(row)
    iRow = row(i);
    iCol = col(i);
    temp = M(iRow,iCol);
    if ~isempty(symvar(temp))
        str1=bracketSubs(char(vpa(temp)));
        str=[var,'(',num2str(iRow),',',num2str(iCol),') = interval(',str1,');'];
    else
        str=[var,'(',num2str(iRow),',',num2str(iCol),') = interval(',char(temp),');']; 
    end
    %write in file
    fprintf(fid, '%s\n', str);
end

%------------- END OF CODE --------------