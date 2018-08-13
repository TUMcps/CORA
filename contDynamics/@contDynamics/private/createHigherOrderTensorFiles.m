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
%    ---
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------   

    % construct auxiliary variables
    N = options.tensorOrder;
    z = [vars.x;vars.u];
    dz = [varsDer.x;varsDer.u];

    % generate the fourth order tensor
    tensor = generateNthTensor(fdyn,z,4);
    func = evalNthTensor(tensor,dz,4);
    func = simplification(func,options);
    str = sprintf('tensor%i_%s',4,name);
    pathFile = [path, filesep, str];
    matlabFunction(func,'File',pathFile,'Vars',{vars.x,vars.u,varsDer.x,varsDer.u});
    disp('tensor 4th-order');
    
    % generate all remaining higher order tensors
    for i = 5:N
        tensor = generateNthTensor(fdyn,z,i,tensor);
        func = evalNthTensor(tensor,dz,i);
        func = simplification(func,options);
        str = sprintf('tensor%i_%s',i,name);
        pathFile = [path, filesep, str];
        matlabFunction(func,'File',pathFile,'Vars',{vars.x,vars.u,varsDer.x,varsDer.u});
        disp(['tensor ',num2str(i),'th-order']);
    end
end

function func = simplification(func,options)
% simplifies the symbolic expression "func" with the specified method
    if isfield(options,'simplify')
        if strcmp(options.simplify,'simplify')
            func  = simplify(func);
        elseif strcmp(options.simplify,'collect')
            func = collect(func,dz);
        elseif ~strcmp(options.simplify,'none')
            error('Wrong value for options.simplify!. Only values ''simplify'', ''collect'' and ''none'' are valid!')
        end
    end
end

%------------- END OF CODE --------------