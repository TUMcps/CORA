function [requiredFiles, requiredFiles_out, allnew, storedData, filepathOld] = ...
    checkTensorRecomputation(obj,fdyn,fcon,fout,vars,options)
% checkTensorRecomputation - checks whether symbolic computations have to
%    be performed or whether the old derivations have remained unchanged
%
% Syntax:
%    [generateFiles, requiredFiles_out, allnew, storedData, filepathOld] = ...
%       checkTensorRecomputation(obj,fdyn,fcon,fout,vars,options)
%
% Inputs:
%    obj - contDynamics object
%    fdyn - symbolic differential equation
%    fcon - symbolic constraint equation (only nonlinDASys)
%    fout - symbolic output equation
%    vars - symbolic variables
%    options - options for reachability analysis, mainly
%                options.tensorOrder and options.lagrangeRem
%
% Outputs:
%    requiredFiles - logical array which tensor files need to be recomputed
%                    (for the dynamic/constaint equations)
%    requiredFiles_out - logical array which tensor files need to be
%                        recomputed (for the output equation)
%    allnew - flag whether all files have to be recomputed
%    storedData - information about computed tensor files
%    filepathOld - path of folder containing tensor files
%
% References:
%    -

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       ---
% Last update:   01-February-2021 (MW, introduction of requiredFiles)
%                18-November-2022 (MW, integrate output equation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set standard path
path = [CORAROOT filesep 'models' filesep 'auxiliary' filesep obj.name];

% init storedData
storedData = [];

% initial assumption: state which files are required
requiredFiles.higherOrders = false(10,1);
% without interval arithmetic
requiredFiles.standard = false(3,1);
% with interval arithmetic
requiredFiles.int = false(3,1);

% copy for output equation...
requiredFiles_out = requiredFiles;

% assume that all files have to be newly computed
allnew = true;

% check required files for dynamic/constaint equation
requiredFiles = aux_checkRequiredFiles(obj,options,requiredFiles,'dynamic');
% check required files for output equation
if all(obj.out_isLinear)
    % only Jacobian required
    requiredFiles_out.standard(1) = true;
else
    requiredFiles_out = aux_checkRequiredFiles(obj,options,requiredFiles_out,'output');
end

% perform a number of checks to see if all / some files need recomputation

% check if dynamics files already exists
filepathOld = [path filesep obj.name '_lastVersion.mat'];
% 1. check if files path exists (derivatives previously computed)
if ~exist(filepathOld,'file')
    return;
end

% load stored data (dynamics and settings) 
load(filepathOld);

% 2. check if something has been previously computed
if ~exist('storedData','var')
    return;
end

% compare stored and actual dynamics (given in symbolic variables)
equalDynEquations = isequal(fdyn,storedData.fdyn);
equalConEquations = isequal(fcon,storedData.fcon);
if isfield(storedData,'fout')
    equalOutEquations = isequal(fout,storedData.fout);
else
    % ...old struct in models/auxiliary
    equalOutEquations = false;
end

% 3. check if equations match
if ~all([equalDynEquations,equalConEquations,equalOutEquations])
    return;
end

% 4. check options.paramInt
if isfield(storedData,'paramInt') ~= isfield(options,'paramInt') ...
    || (isfield(storedData,'paramInt') && (isa(storedData.paramInt,'interval') ~= isa(options.paramInt,'interval'))) ...
    || (isfield(storedData,'paramInt') && isa(storedData.paramInt,'interval') && ~isequal(storedData.paramInt,options.paramInt))
    % existence of paramInt, class and set (interval) have to be equal
        return;
end

% check options.lagrangeRem: assign default if non-existant
if isfield(options,'lagrangeRem')
    temp1 = options.lagrangeRem;
else
    temp1.simplify = getDefaultValue('lagrangeRem.simplify',obj,struct(),options,'options');
    temp1.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',obj,struct(),options,'options');
    temp1.method = getDefaultValue('lagrangeRem.method',obj,struct(),options,'options');
end
if isfield(storedData,'lagrangeRem')
    temp2 = storedData.lagrangeRem;
else
    temp2.simplify = getDefaultValue('lagrangeRem.simplify',obj,struct(),options,'options');
    temp2.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',obj,struct(),options,'options');
    temp2.method = getDefaultValue('lagrangeRem.method',obj,struct(),options,'options');
end

% 5. check setting 'simplify':
if ~isfield(temp1,'simplify')
    temp1.simplify = getDefaultValue('lagrangeRem.simplify',obj,struct(),options,'options');
end
if ~isfield(temp2,'simplify')
    temp2.simplify = getDefaultValue('lagrangeRem.simplify',obj,struct(),options,'options');
end
if ~strcmp(temp1.simplify,temp2.simplify)
    return;
end

% 6. check setting 'tensorParallel':
if (isfield(temp1,'tensorParallel') && isfield(temp2,'tensorParallel') ...
        && temp1.tensorParallel ~= temp2.tensorParallel) ...
        || (isfield(temp1,'tensorParallel') && temp1.tensorParallel) ...
        || (isfield(temp2,'tensorParallel') && temp2.tensorParallel)
    return;
end

% 7. check setting 'method':
if ~isfield(temp1,'method')
    temp1.method = getDefaultValue('lagrangeRem.method',obj,struct(),options,'options');
end
if ~isfield(temp2,'method')
    temp2.method = getDefaultValue('lagrangeRem.method',obj,struct(),options,'options');
end
if ~strcmp(temp1.method,temp2.method)
    return;
end

% 8. check setting 'replacements'
if isfield(temp1,'replacements') ~= isfield(temp2,'replacements') ...
        || (isfield(temp1,'replacements') && ~isempty(vars.p) ...
        && ~isequal(temp1.replacements(vars.x,vars.u,vars.p),temp2.replacements)) ...
        || (isfield(temp1,'replacements') && isempty(vars.p) ...
        && ~isequal(temp1.replacements(vars.x,vars.u),temp2.replacements))
    return;
end

% here: all checks reveal equality between previous and current setting
%       -> remove already computed tensor files from required ones
allnew = false;

% names of tensor files
prefixes = {'jacobian','hessianTensor','thirdOrderTensor'};

% dynamic equation
for i=1:3
    % not required anymore if already computed
    if isfile([path filesep prefixes{i} '_' obj.name '.m'])
        requiredFiles.standard(i) = false;
    end
    if isfile([path filesep prefixes{i} 'Int_' obj.name '.m'])
        requiredFiles.int(i) = false;
    end
end
% no separation of file names (with Int) for 4th or higher orders
% for i=4:10
%     prefixes{i} = ['tensor' num2str(i)];
%     if isfile([path filesep prefixes{i} '_' obj.name '.m'])
%         requiredFiles.higherOrders(i) = false;
%     end
% end

% output equation
for i=1:3
    % not required anymore if already computed
    if isfile([path filesep 'out_' prefixes{i} '_' obj.name '.m'])
        requiredFiles_out.standard(i) = false;
    end
    if isfile([path filesep 'out_' prefixes{i} 'Int_' obj.name '.m'])
        requiredFiles_out.int(i) = false;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function requiredFiles = aux_checkRequiredFiles(obj,options,requiredFiles,eq)

if ( isa(obj,'nonlinearSys') || isa(obj,'nonlinDASys') || isa(obj,'nonlinearSysDT') ) ...
        && isfield(options,'alg') && contains(options.alg,'adaptive')
    % adaptive algorithm (only nonlinearSys: options.alg = 'lin|poly-adaptive')
    if contains(options.alg,'lin')
        requiredFiles.standard(1:2) = true;
        requiredFiles.int(2:3) = true;
    else
        requiredFiles.standard(1:2) = true;
        requiredFiles.int(3) = true;
    end
else
    if strcmp(eq,'dynamic')
        tensorOrder = options.tensorOrder;
    elseif strcmp(eq,'output')
        tensorOrder = options.tensorOrderOutput;
    end
    if tensorOrder == 1
        requiredFiles.standard(1) = true;
    elseif tensorOrder <= 3
        requiredFiles.standard(1:tensorOrder-1) = true;
        requiredFiles.int(tensorOrder) = true;
    else
        requiredFiles.standard(1:3) = true;
        requiredFiles.higherOrders(4:tensorOrder) = true;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
