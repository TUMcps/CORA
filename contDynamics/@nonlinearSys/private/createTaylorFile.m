function createTaylorFile(fTaylor,path,name)
% createTaylorFile - generates an mFile of the Taylor expansion of the
% nonlinear systems object
%
% Syntax:  
%    createTaylorFile(obj,fTaylor,path,name)
%
% Inputs:
%    obj - nonlinear sytsem object
%    fTaylor - Taylor expressions
%    path - path for file location
%    name - name of the corresponding nonlinear system
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
% Written:      06-December-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


fid = fopen([path '/taylorModel_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function f=taylorModel_',name,'(x,u)']);
for k=1:length(fTaylor)
    str=['f(',num2str(k),',1)=',char(fTaylor(k,1)),';'];
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
        
    %write in file
    fprintf(fid, '%s\n', str);
end

%close file
fclose(fid);

%------------- END OF CODE --------------