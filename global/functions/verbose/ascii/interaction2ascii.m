function interaction2ascii
% interaction2ascii - ???
%
% Syntax:
%    interaction2ascii
%
% Inputs:
%    - 
%
% Outputs:
%    - 
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       30-May-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%open .mat file
[FileName,PathName] = uigetfile();
cd(PathName);
file=load(FileName);
Omega=file.Omega;

[nrOfAmodes, nrOfBmodes]=size(Omega);

for carAmode=1:nrOfAmodes
    for carBmode=1:nrOfBmodes
        %create filename
        fName=strrep(FileName,'.mat',[num2str(carAmode),...
            '_',num2str(carBmode),'.txt']);
        %open file for writing
        fid = fopen(fName,'w');
        %write to file
        aux_saveMatrix(fid,Omega{carAmode,carBmode});
        %close file
        status = fclose(fid)
    end
end

end


% Auxiliary functions -----------------------------------------------------

function aux_saveMatrix(fid,Omega)
%1. value: nr of rows
%2. value: nr of columns
%next values: elements values of the matrix

    %get nr of rows and columns
    [nrOfRows,nrOfColumns]=size(Omega);
    
    %save nr of rows, columns
    fprintf(fid, '%3i %3i \n', [nrOfRows, nrOfColumns]);
    
    %save nonzero elements: row, coulmn, value
    for iRow=1:nrOfRows
        for iCol=1:nrOfColumns
            fprintf(fid, '%0.4f \n', Omega(iRow,iCol));
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
