function createParametricDynamicFile(obj,mFile,path,name)
% createParametricDynamicFile - generates an mFile of the dynamic equations
%    sorted by parameter influences
%
% Syntax:
%    createParametricDynamicFile(obj,mFile,path,name)
%
% Inputs:
%    obj - nonlinear system object
%    mFile - function handle for the function (dynamic or output)
%    path - file-path to the folder containing the model files
%    name - function name for the parametric dynamic file
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
% Written:       01-June-2011
% Last update:   02-June-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%create symbolic variables
vars = symVariables(obj,'LRbrackets');

%insert symbolic variables into the system equations
f=mFile(vars.x,vars.u,vars.p);

%init
fcell=cell(1,obj.nrOfParam+1);
%part without parameters
fcell{1} = subs(f,vars.p,zeros(obj.nrOfParam,1));
%part with parameters
I = eye(obj.nrOfParam); %identity matrix
for i=1:obj.nrOfParam
    fcell{i+1} = subs(f,vars.p,I(:,i)) - fcell{1};
end


fid = fopen([path filesep name '.m'],'w');
fprintf(fid, '%s\n\n', ['function f=' name '(x,u)']);
for k=1:length(fcell)
    for i=1:length(f)
        str=['f{',num2str(k),'}(',num2str(i),',1)=',char(fcell{k}(i,1)),';'];
        str=bracketSubs(str);
        %write in file
        fprintf(fid, '%s\n', str);
    end
end

%close file
fclose(fid);

% ------------------------------ END OF CODE ------------------------------
