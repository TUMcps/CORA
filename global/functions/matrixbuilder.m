function [M]=matrixbuilder(rows,columns,type)
% updated: 18-June-2009, MA
% updated: 06-November-2009, MA

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
        M(iColumn,rows*(iColumn-1)+iRow+1)=1;
    end
end
