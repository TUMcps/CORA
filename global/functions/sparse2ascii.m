function sparse2ascii
%last update: 28.08.08, MA

%create dummy Markov chain to avoid errors while loading
MCdummy=markovchain([],[]);
fieldDummy=partition([1,2],[2]);

%open .mat file
[FileName,PathName] = uigetfile();
cd(PathName);
file=load(FileName);
MC=file.MC;

T=get(MC,'T');

nrOfModes=length(T.T);

for carMode=1:nrOfModes
    %time point solution
    %create filename
    fName=strrep(FileName,'.mat',['_T_',num2str(carMode),'.txt']);
    %open file for writing
    fid = fopen(fName,'w');
    %write to file
    saveMatrix(fid,T.T{carMode});
    %close file
    status = fclose(fid)
    
    %time interval solution
    %create filename
    fName=strrep(FileName,'.mat',['_OT_',num2str(carMode),'.txt']);
    %open file for writing
    fid = fopen(fName,'w');
    %write to file
    saveMatrix(fid,T.OT{carMode});
    %close file
    status = fclose(fid)    
end


%--------------------------------------------------------------------------
function saveMatrix(fid,T)
%1. value: nr of rows
%2. value: nr of columns
%3. value: nr of non-zero elements
%next values: row, column, value

%get nr of rows and columns
[nrOfRows,nrOfColumns]=size(T);

%get non-zero rows and columns
[nonZeroRow,nonZeroCol]=find(T);

%get nr of non-zero elements
nrOfNonZeroElements=length(nonZeroRow);

%save nr of rows, columns, non-zero elements
fprintf(fid, '%3i %3i %3i \n', [nrOfRows, nrOfColumns, nrOfNonZeroElements]);

%save nonzero elements: row, column, value
for i=1:nrOfNonZeroElements 
    val=T(nonZeroRow(i),nonZeroCol(i));
    fprintf(fid, '%3i %3i %0.4f \n', [nonZeroRow(i),nonZeroCol(i), val]);
end