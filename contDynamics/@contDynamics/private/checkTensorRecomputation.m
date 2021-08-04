function [requiredFiles, storedData, filepathOld] = checkTensorRecomputation(obj,fdyn,fcon,vars,options)
% checkTensorRecomputation - checks whether symbolic computations have to
%    be performed or whether the old derivations have remained unchanged
%
% Syntax:  
%    [generateFiles, storedData, filepathOld] = checkTensorRecomputation(obj,fdyn,fcon,vars,options)
%
% Inputs:
%    obj - contDynamics object
%    fdyn - symbolic differential equation
%    fcon - symbolic constraint equation (only nonlinDASys)
%    vars - symbolic variables
%    options - options for reachability analysis, mainly
%                options.tensorOrder and options.lagrangeRem
%
% Outputs:
%    requiredFiles - logical array which tensor files need to be recomputed
%    storedData - information about computed tensor files
%    filepathOld - path of folder containing tensor files
%
% References:
%    -

% Author:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:      ---
% Last update:  01-February-2021 (MW, introduction of requiredFiles)
% Last revision:---

%------------- BEGIN CODE --------------

% set standard path
path = [coraroot filesep 'models' filesep 'auxiliary' filesep obj.name];

% init storedData
storedData = [];

% initial assumption: state which files are required
requiredFiles.higherOrders = false(10,1);
requiredFiles.standard = false(3,1);    % without interval arithmetic
requiredFiles.int = false(3,1);         % with interval arithmetic
if isa(obj,'nonlinearSys') && contains(options.alg,'adaptive')
    % adaptive algorithm (only nonlinearSys: options.alg = 'lin|poly-adaptive')
    if contains(options.alg,'lin')
        requiredFiles.standard(1:2) = true;
        requiredFiles.int(2:3) = true;
    else
        requiredFiles.standard(1:2) = true;
        requiredFiles.int(3) = true;
    end
else
    if options.tensorOrder <= 3
        requiredFiles.standard(1:options.tensorOrder-1) = true;
        requiredFiles.int(options.tensorOrder) = true;
    else
        requiredFiles.standard(1:3) = true;
        requiredFiles.higherOrders(4:options.tensorOrder) = true;
    end
end
requiredFiles.allnew = true;

% perform a number of checks to see if all / some files need recomputation

% check if dynamics files already exists
filepathOld = [path filesep obj.name '_lastVersion.mat'];
% 1. check if files path exists (derivatives previously computed)
if ~exist(filepathOld,'file')
    return;
end

% load stored data (dynamics and settings) 
load(filepathOld);

% 2. compare stored and actual dynamics (given in symbolic variables)
if ~exist('storedData','var')
    return;
end

% compare stored and actual dynamics (given in symbolic variables)
equalDynEquations = isequal(fdyn,storedData.fdyn);
equalConEquations = isequal(fcon,storedData.fcon);

% 3. check if equations match
if ~all([equalDynEquations,equalConEquations])
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
    temp1.simplify = getDefaultValue('lagrangeRem.simplify',obj,'options');
    temp1.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',obj,'options');
    temp1.method = getDefaultValue('lagrangeRem.method',obj,'options');
end
if isfield(storedData,'lagrangeRem')
    temp2 = storedData.lagrangeRem;
else
    temp2.simplify = getDefaultValue('lagrangeRem.simplify',obj,'options');
    temp2.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',obj,'options');
    temp2.method = getDefaultValue('lagrangeRem.method',obj,'options');
end

% 5. check setting 'simplify':
if ~isfield(temp1,'simplify')
    temp1.simplify = getDefaultValue('lagrangeRem.simplify',obj,'options');
end
if ~isfield(temp2,'simplify')
    temp2.simplify = getDefaultValue('lagrangeRem.simplify',obj,'options');
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
    temp1.method = getDefaultValue('lagrangeRem.method',obj,'options');
end
if ~isfield(temp2,'method')
    temp2.method = getDefaultValue('lagrangeRem.method',obj,'options');
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
requiredFiles.allnew = false;

% names of tensor files
prefixes = {'jacobian','hessianTensor','thirdOrderTensor'};
% not required anymore if already computed
for i=1:3
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

end

%------------- END OF CODE --------------
