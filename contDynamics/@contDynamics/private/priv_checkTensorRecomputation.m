function [requiredFiles,requiredFiles_out,requiredData,deleteAll] = ...
    priv_checkTensorRecomputation(sys,fdyn,fcon,fout,path,options)
% priv_checkTensorRecomputation - checks whether symbolic computations have to
%    be performed or whether the old derivations have remained unchanged
%
% Syntax:
%    [generateFiles,requiredFiles_out,requiredData] = ...
%       priv_checkTensorRecomputation(sys,fdyn,fcon,fout,path,options)
%
% Inputs:
%    sys - contDynamics object
%    fdyn - symbolic differential equation
%    fcon - symbolic constraint equation (only nonlinDASys)
%    fout - symbolic output equation
%    path - path to folder where tensor files are computed
%    options - options for reachability analysis
%           .tensorOrder
%           .lagrangeRem
%
% Outputs:
%    requiredFiles - logical array which tensor files need to be recomputed
%                    (for the dynamic/constaint equations)
%    requiredFiles_out - logical array which tensor files need to be
%                        recomputed (for the output equation)
%    requiredData - information about settings for tensor file generation
%    deleteAll - true/false whether all files should be deleted
%
% Example:
%    -

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       ---
% Last update:   01-February-2021 (MW, introduction of requiredFiles)
%                18-November-2022 (MW, integrate output equation)
% Last revision: 07-October-2024 (MW, complete refactor)

% ------------------------------ BEGIN CODE -------------------------------

% list required files for dynamic/constraint equation and output equation
% given the current dynamics and tensorOrder
requiredFiles = aux_listRequiredFiles(sys,true,options);
requiredFiles_out = aux_listRequiredFiles(sys,false,options);
deleteAll = false;

% check if expected path to .mat-file containing information about the 
% folder content exists (true if derivatives previously computed); if yes,
% try loading stored data from that file
try
    load([path filesep sys.name '_lastVersion.mat'],'storedData');
catch ME
    storedData = struct([]);
end

% generate equivalent struct to storedData with current settings
requiredData = aux_requiredData(sys,fdyn,fcon,fout,options);

% also save current CORAversion and datetime of generation
requiredData.CORAversion = CORAVERSION;
requiredData.timeStamp = datetime;

% structs must be equal, otherwise we must not check for files that are
% computed already (see below) since we have different dynamics/settings!
if isempty(storedData) || ~aux_compareData(storedData,requiredData)
    deleteAll = true;
    return;
end

% update required files: remove those that are computed already
[requiredFiles,requiredFiles_out] = aux_updateRequiredFiles(sys,path,...
    requiredFiles,requiredFiles_out);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_compareData(S1,S2)
% compares the data between the two structs; normally, isequal(S1,S2)
% should suffice, but .lagrangeRem.replacements is a function handle, for
% which the following behavior holds true:
%   f = @(x) x(1)*x(2);
%   g = @(x) x(1)*x(2);
%   assert(~isequal(f,g));  % <-- this is a problem for us
% hence, we compare the struct while omitting any occurrences of the field
% .lagrangeRem.replacements and compare them separately

% both must have .lagrangeRem -> check lagrangeRem.replacements separately
if xor(isfield(S1.lagrangeRem,'replacements'),isfield(S2.lagrangeRem,'replacements'))
    res = false;
    return
end

if isfield(S1.lagrangeRem,'replacements') && isfield(S2.lagrangeRem,'replacements')
    if ~isequalFunctionHandle(S1.lagrangeRem.replacements,S2.lagrangeRem.replacements)
        res = false;
        return
    end
    % remove replacements from further checking
    S1.lagrangeRem = rmiffield(S1.lagrangeRem,'replacements');
    S2.lagrangeRem = rmiffield(S2.lagrangeRem,'replacements');
end

% unfortunately, we also need to check paramInt separately... the reason is
% that the isequal-call below for the struct (correctly!) uses the
% isequal-function for interval or double. In this case, a comparison of
% interval(a,a) and a would return true (as intended). However, we must
% ensure that both .paramInt are of the same class since the computation of
% derivatives differs whether .paramInt is an interval or a double.
if xor(isfield(S1.lagrangeRem,'replacements'),isfield(S2.lagrangeRem,'replacements'))
    res = false;
    return
end
if isfield(S1,'paramInt') && isfield(S2,'paramInt') ...
        && xor(isa(S1.paramInt,'interval'),isa(S2.paramInt,'interval'))
    res = false;
    return    
end

% remove version and generated as they do not need to be equal
S1 = rmiffield(S1,{'CORAversion','timeStamp'});
S2 = rmiffield(S2,{'CORAversion','timeStamp'});

% use built-in broadcast of isequal functions via struct comparison
res = isequal(S1,S2);

end

function requiredData = aux_requiredData(sys,fdyn,fcon,fout,options)

% save dynamics/constraint/output equation
requiredData.fdyn = fdyn;
requiredData.fcon = fcon;
requiredData.fout = fout;

% save settings for Lagrange remainder: here, we are a bit lenient with the
% caller... add default values for all non-set fields
if ~isfield(options,'lagrangeRem')
    requiredData.lagrangeRem.simplify = getDefaultValue('lagrangeRem.simplify',sys,struct(),options,'options');
    requiredData.lagrangeRem.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',sys,struct(),options,'options');
    requiredData.lagrangeRem.method = getDefaultValue('lagrangeRem.method',sys,struct(),options,'options');
    % no .replacements here
else
    requiredData.lagrangeRem = options.lagrangeRem;
    if ~isfield(options.lagrangeRem,'simplify')
        requiredData.lagrangeRem.simplify = getDefaultValue('lagrangeRem.simplify',sys,struct(),options,'options');
    end
    if ~isfield(options.lagrangeRem,'tensorParallel')
        requiredData.lagrangeRem.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',sys,struct(),options,'options');
    end
    if ~isfield(options.lagrangeRem,'method')
        requiredData.lagrangeRem.method = getDefaultValue('lagrangeRem.method',sys,struct(),options,'options');
    end
    % no .replacements here...
end

% only for dynamical systems with parameters
if isfield(options,'paramInt')
    requiredData.paramInt = options.paramInt;
end

end

function files = aux_listRequiredFiles(sys,isState,options)
% list the files that need to be generated according to the given settings
%   isState == true    for state equation (differential equation)
%   isState == false   for output equation 

% initial assumption: state which files are required
files.higherOrders = false(10,1);
% without interval arithmetic
files.standard = false(3,1);
% with interval arithmetic
files.int = false(3,1);

% shortcut if we deal with the output equation
if ~isState
    if isa(sys,'nonlinearARX')
        % no output files required
        return
    elseif all(sys.out_isLinear)
        % only Jacobian required... all further derivatives are zero
        files.standard(1) = true;
        return
    end
end

% conformant synthesis, no output equation
if isfield(options,'cs')
    files.standard(1) = true;
    return
end

% adaptive algorithm... requires multiple no-interval-arithmetic and
% interval arithmetic tensors
if ( isa(sys,'nonlinearSys') || isa(sys,'nonlinDASys') || isa(sys,'nonlinearSysDT') ) ...
        && isfield(options,'alg') && contains(options.alg,'adaptive')
    % only nonlinearSys has 'lin' and 'poly' in options.alg
    if contains(options.alg,'lin')
        files.standard(1:2) = true;
        files.int(2:3) = true;
    else  % 'poly'
        files.standard(1:2) = true;
        files.int(3) = true;
    end
    return
end

% read out tensor order
if isState
    tensorOrder = options.tensorOrder;
else
    tensorOrder = options.tensorOrderOutput;
end

% always compute Jacobian...
files.standard(1) = true;
% standard requirements for derivatives computation
if tensorOrder <= 3
    files.standard(1:tensorOrder-1) = true;
    files.int(tensorOrder) = true;
else
    files.standard(1:3) = true;
    files.higherOrders(4:tensorOrder) = true;
end

end

function [requiredFiles,requiredFiles_out] = aux_updateRequiredFiles(...
    sys,path,requiredFiles,requiredFiles_out)

% (hardcoded) names of tensor files
prefixes = {'jacobian','hessianTensor','thirdOrderTensor'};

% go through files in path and find which ones are there already...
for i=1:3
    % dynamic equation
    requiredFiles.standard(i) = requiredFiles.standard(i) ...
        && ~isfile([path filesep prefixes{i} '_' sys.name '.m']);
    requiredFiles.int(i) = requiredFiles.int(i) ...
        && ~isfile([path filesep prefixes{i} 'Int_' sys.name '.m']);
    % output equation
    requiredFiles_out.standard(i) = requiredFiles_out.standard(i) ...
        && ~isfile([path filesep 'out_' prefixes{i} '_' sys.name '.m']);
    requiredFiles_out.int(i) = requiredFiles_out.int(i) ...
        && ~isfile([path filesep 'out_' prefixes{i} 'Int_' sys.name '.m']);
end

% no separation of file names (with Int) for 4th or higher orders
for i=4:10
    prefixes{i} = ['tensor' num2str(i)];
    requiredFiles.higherOrders(i) = requiredFiles.higherOrders(i) ...
        && ~isfile([path filesep prefixes{i} '_' sys.name '.m']);
end

end

% ------------------------------ END OF CODE ------------------------------
