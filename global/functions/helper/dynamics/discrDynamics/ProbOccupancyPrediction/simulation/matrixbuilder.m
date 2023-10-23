function M = matrixbuilder(rows,columns,type)
% matrixbuilder - ???
%
% Syntax:
%    M = matrixbuilder(rows,columns,type)
%
% Inputs:
%    rows - ???
%    columns - ???
%    type - ???
%
% Outputs:
%    M - ???
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       18-June-2009
% Last update:   06-November-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if type==0
    %init M
    M=sparse(rows,rows*columns+1);
    for iRow=1:rows
        iColumn=1:columns;
        M(iRow,iRow+rows*(iColumn-1)+1)=1;
    end
else
    %init M
    M=sparse(columns,rows*columns+1);    
    for iColumn=1:columns
        iRow=1:rows;
        M(iColumn,iRow+rows*(iColumn-1)+1)=1;
    end
end

% ------------------------------ END OF CODE ------------------------------
