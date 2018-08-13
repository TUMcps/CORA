function writeSparseMatrix(M,var,fid)


%write each row
[row,col] = find(M~=0);

for i=1:length(row)
    iRow = row(i);
    iCol = col(i);
    str=bracketSubs(char(M(iRow,iCol)));
    str=[var,'(',num2str(iRow),',',num2str(iCol),') = ',str,';'];
    %write in file
    fprintf(fid, '%s\n', str);
end