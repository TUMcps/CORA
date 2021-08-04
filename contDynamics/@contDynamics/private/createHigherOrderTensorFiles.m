function createHigherOrderTensorFiles(fdyn,vars,varsDer,path,name,options)
% createHigherOrderTensorFiles - create tensor files with order > 3
%
% Syntax:  
%    createHigherOrderTensorFiles(fdyn,vars,varsDer,path,name,options)
%
% Inputs:
%    fdyn - symbolic function
%    vars - struct containing the symbolic variables of the function
%    varsDer - struct containing the symbolic derivatives of the variables
%    path - file-path to the folder where the generated files are stored
%    name - name of the dynamical system
%    options - struct containing the algorithm options
%
% Outputs:
%    -
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Niklas Kochdumper
% Written:      08-February-2018
% Last update:  02-February-2021 (MW, remove code duplicates)
% Last revision:---

%------------- BEGIN CODE --------------   

    % construct auxiliary variables
    N = options.tensorOrder;
    z = [vars.x;vars.u];
    dz = [varsDer.x;varsDer.u];
    
    tensor = []; % init for first call of generateNthTensor (order 4)
    % generate all higher-order tensors
    for i = 4:N
        tensor = generateNthTensor(fdyn,z,i,tensor);
        func = evalNthTensor(tensor,dz,i);
        func = simplification(func,options,dz);
        str = sprintf('tensor%i_%s',i,name);
        pathFile = [path, filesep, str];
        if ~isempty(vars.p)
            matlabFunction(func,'File',pathFile,'Vars',{vars.x,vars.u,varsDer.x,varsDer.u,vars.p});
        else
            matlabFunction(func,'File',pathFile,'Vars',{vars.x,vars.u,varsDer.x,varsDer.u});
        end
        disp(['... compute symbolic tensor ' num2str(i) 'th-order']);
    end
    
end


function func = simplification(func,options,dz)
% simplifies the symbolic expression "func" with the specified method
    
    if isfield(options,'lagrangeRem')
        temp = options.lagrangeRem;
        if isfield(temp,'simplify')
            if strcmp(temp.simplify,'simplify')
                func  = simplify(func);
            elseif strcmp(temp.simplify,'collect')
                func = collect(func,dz);
            end
        end
    end
    
end

%------------- END OF CODE --------------