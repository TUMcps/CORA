function empty = writeSparseMatrix(M,var,fid)
% writeSparseMatrix - write a sparse matrix to file
%
% Syntax:  
%    empty = writeSparseMatrix(M,var,fid)
%
% Inputs:
%    M - symbolic matrix
%    var - name of the matrix that is written
%    fid - identifier of the file to which the matrix is written
%
% Outputs:
%    empty - 1 if matrix is empty, 0 if not
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

% write each row
[row,col] = find(M~=0);

if ~isempty(row)

    empty = 0;
    
    for i=1:length(row)
        iRow = row(i);
        iCol = col(i);
        str=bracketSubs(char(M(iRow,iCol)));
        str=[var,'(',num2str(iRow),',',num2str(iCol),') = ',str,';'];
        % write in file
        fprintf(fid, '%s\n', str);
    end
    
else
    empty = 1;
end

%------------- END OF CODE --------------