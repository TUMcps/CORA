function createFileFromFunction(f,name,output,input)
% createFileFromFunction - generates an m-file of a specified function
%
% Syntax:
%    createFileFromFunction(f,name,output,input)
%
% Inputs:
%    f - function
%    name - name of function
%    output - symbol for the output
%    input - string of comma-separated inputs 
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
% See also: ---

% Authors:       Matthias Althoff
% Written:       16-April-2022
% Last update:   25-April-2023 (MW, update path)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];
if ~exist(path,'dir')
   mkdir(path); 
end
addpath(path);

% create file
fid = fopen([path filesep name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function ',output,' = ',name,'(',input,')']);

if ~isempty(f)
    % non-empty f
    for k=1:length(f)
        str=[output,'(',num2str(k),',1)=',char(f(k,1)),';'];
        %generate left and right brackets
        str=strrep(str,'L','(');
        str=strrep(str,'R',')');

        %write in file
        fprintf(fid, '%s\n', str);
    end
else
    % empty f
    fprintf(fid, '%s\n', 'f=[];');
end

%close file
fclose(fid);

% remove and add path so that file can be found
warOrig = warning;
warning('off','all');
rmpath(path);
warning(warOrig);
addpath(path);

% ------------------------------ END OF CODE ------------------------------
