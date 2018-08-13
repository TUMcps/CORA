function matrix2ascii(matrix)
%last update: 29.10.08, MA


%create filename
fName=[pwd,'/matrix.txt'];
%open file for writing
fid = fopen(fName,'w');

%write to file
%get nr of rows and columns
[nrOfRows,nrOfColumns]=size(matrix);

%save nr of rows, columns
fprintf(fid, '%3i %3i \n', [nrOfRows, nrOfColumns]);

%save nonzero elements: row, coulmn, value
for iRow=1:nrOfRows
    for iCol=1:nrOfColumns
        fprintf(fid, '%0.4f \n', matrix(iRow,iCol));
    end
end

%close file
status = fclose(fid);