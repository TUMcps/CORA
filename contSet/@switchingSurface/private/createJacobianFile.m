function createJacobianFile(J,path,name)
% createJacobianFile - generates an mFile that allows to compute the
% jacobian of the intersectionTimeModel
%
% Syntax:  
%    createJacobianFile(J,path)
%
% Inputs:
%    J - symbolic jacobian
%    path - file path
%    name - file name extension
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

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%set directory
cd(path) %<-- change directory here

fid = fopen(['jacobian_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function [J]=jacobian_',name,'(A,u,x_c,x_0,y)']);

% write "J=["
fprintf(fid, '%s', 'J=[');
% write rest of matrix
writeMatrix(J,fid);


%close file
fclose(fid);


function writeMatrix(M,fid)

%write each row
for iRow=1:(length(M(:,1)))
    if (length(M(1,:))-1)>0
        for iCol=1:(length(M(1,:))-1)
            str=charSubs(char(M(iRow,iCol)));
            str=[str,','];
            %write in file
            fprintf(fid, '%s', str);
        end
    else
        iCol = 0; %for vectors
    end
    if iRow<length(M(:,1))
        %write last element
        str=charSubs(char(M(iRow,iCol+1)));
        str=[str,';...'];
        %write in file
        fprintf(fid, '%s\n', str);
    else
        %write last element
        str=charSubs(char(M(iRow,iCol+1)));
        str=[str,'];'];
        %write in file
        fprintf(fid, '%s\n\n', str);   
    end
end


function [str]=charSubs(str)

%generate left and right brackets
str=strrep(str,'L','(');
str=strrep(str,'R',')');
str=strrep(str,'CO',',');


%------------- END OF CODE --------------