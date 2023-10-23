function writeMatrix(M,fid)
% writeMatrix - ???
%
% Syntax:
%    writeMatrix(M,fid)
%
% Inputs:
%    M - symbolic matrix
%    fid - identifier of the file to which the matrix is written
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

%write each row
for iRow=1:(length(M(:,1)))
    if (length(M(1,:))-1)>0
        for iCol=1:(length(M(1,:))-1)
            str=bracketSubs(char(M(iRow,iCol)));
            str=[str,','];
            %write in file
            fprintf(fid, '%s', str);
        end
    else
        iCol = 0; %for vectors
    end
    if iRow<length(M(:,1))
        %write last element
        str=bracketSubs(char(M(iRow,iCol+1)));
        str=[str,';...'];
        %write in file
        fprintf(fid, '%s\n', str);
    else
        %write last element
        str=bracketSubs(char(M(iRow,iCol+1)));
        str=[str,'];'];
        %write in file
        fprintf(fid, '%s\n\n', str);   
    end
end

% ------------------------------ END OF CODE ------------------------------
