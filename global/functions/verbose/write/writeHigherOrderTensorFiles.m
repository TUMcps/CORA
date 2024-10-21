function writeHigherOrderTensorFiles(fdyn,vars,varsDer,path,fname,options)
% writeHigherOrderTensorFiles - create tensor files with order > 3
%
% Syntax:
%    writeHigherOrderTensorFiles(fdyn,vars,varsDer,path,fname,options)
%
% Inputs:
%    fdyn - symbolic function
%    vars - struct containing the symbolic variables of the function
%    varsDer - struct containing the symbolic derivatives of the variables
%    path - file-path to the folder where the generated files are stored
%    fname - name of the dynamical system
%    options - struct containing the algorithm options
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
% See also: none

% Authors:       Niklas Kochdumper
% Written:       08-February-2018
% Last update:   02-February-2021 (MW, remove code duplicates)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% concatenate state and inputs
z = [vars.x;vars.u];
dz = [varsDer.x;varsDer.u];

% init for first call of generateNthTensor (order 4)
tensor = [];

% generate all higher-order tensors
for i = 4:options.tensorOrder
    tensor = generateNthTensor(fdyn,z,i,tensor);
    func = evalNthTensor(tensor,dz,i);
    func = aux_simplification(func,options,dz);
    str = sprintf('tensor%i_%s',i,fname);
    pathFile = [path, filesep, str];
    matlabFunction(func,'File',pathFile,'Vars',...
        {vars.x,vars.u,varsDer.x,varsDer.u,vars.p});
    if options.verbose
        disp(['... compute symbolic tensor ' num2str(i) 'th-order']);
    end
end
    
end


% Auxiliary functions -----------------------------------------------------

function func = aux_simplification(func,options,dz)
% simplifies the symbolic expression "func" with the specified method
    
    if isfield(options,'lagrangeRem')
        if isfield(options.lagrangeRem,'simplify')
            if strcmp(options.lagrangeRem.simplify,'simplify')
                func = simplify(func);
            elseif strcmp(options.lagrangeRem.simplify,'collect')
                func = collect(func,dz);
            end
        end
    end
    
end

% ------------------------------ END OF CODE ------------------------------
