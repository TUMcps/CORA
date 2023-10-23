function createModelFactors(c)
% createModelFactors - creates file returning model factors
%
% Syntax:
%    createModelFactors(c)
%
% Inputs:
%    c - symbolic expressions of model factors
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       04-May-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

path = [CORAROOT filesep 'models' filesep 'auxiliary'];
cd(path);

fid = fopen('modelFactors.m','w');
fprintf(fid, '%s\n\n', 'function c=modelFactors(u)');
for k=1:length(c(:,1))
    str=['c(',num2str(k),',1)=',char(c(k,1)),';'];
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
    
    %write in file
    fprintf(fid, '%s\n', str);
end

%close file
status = fclose(fid);

% ------------------------------ END OF CODE ------------------------------
