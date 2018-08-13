function createIntersectionTimeModelFile(eq,path,name)
% createIntersectionTimeModelFile - generates an mFile that represents the 
% intersection time model for a given switching surface h(x) = 0
%
% Syntax:  
%    createIntersectionTimeModelFile(eq,path,name)
%
% Inputs:
%    eq - symbolic equation
%    path - path of file
%    name - file name 
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

fid = fopen(['intersectionTimeModel_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function h = intersectionTimeModel_',name,'(A,u,x_c,x_0,y)']);

if ~isempty(eq)
    str=['h=',char(eq),';'];
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
    str=strrep(str,'CO',',');

    %write in file
    fprintf(fid, '%s\n', str);
else
    %write in file
    fprintf(fid, '%s\n', 'h=[];');
end

%close file
fclose(fid);


%------------- END OF CODE --------------