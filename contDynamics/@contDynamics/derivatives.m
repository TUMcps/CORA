function derivatives(sys,varargin)
% derivatives - computes multivariate derivatives (Jacobians, Hessians, 
%    third-order tensors, higher-order tensors) of nonlinear systems
%    symbolically; the result is stored in m-files to which the
%    contDynamics object has access via handles in its properties
%
% Syntax:
%    derivatives(sys)
%    derivatives(sys,options)
%    derivatives(sys,options,path)
%
% Inputs:
%    sys - contDynamics object
%    options - options struct, with fields
%       .tensorOrder
%       .tensorOrderOutput
%       .lagrangeRem.simplify
%       .lagrangeRem.tensorParallel
%       .lagrangeRem.method
%       .lagrangeRem.replacements
%       .verbose: output text on command window
%    path - filepath to save generated files
%
% Outputs:
%    -
%
% Example:
%    f = @(x,u) [x(1) - u(1); x(1)*x(2)];
%    sys = nonlinearSys('nln',f);
%    derivatives(sys);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       29-October-2007 (MA)
% Last update:   07-September-2012 (MA)
%                12-October-2015 (MA)
%                08-April-2016 (MA)
%                10-June-2017 (NK)
%                15-June-2017 (NK)
%                28-June-2017 (NK)
%                06-July-2017
%                15-July-2017 (NK)
%                16-July-2017 (NK)
%                05-November-2017 (MA, generalization for all contDynamics classes)
%                12-November-2017 (MA)
%                03-December-2017 (MA)
%                14-January-2018 (MA)
%                12-November-2018 (NK, removed lagrange remainder files)
%                29-January-2021 (MW, simplify checks, outsource tensor check)
%                18-November-2022 (MW, integrate output equation)
%                22-March-2024 (LL, add verbose option)
% Last revision: 07-October-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,3);
defaultOptions = aux_defaultOptions(sys);
defaultPath = [CORAROOT filesep 'models' filesep 'auxiliary' filesep sys.name];
[options,path] = setDefaultValues({defaultOptions,defaultPath},varargin);

% be lenient with verbose output...
options.verbose = isfield(options,'verbose') && options.verbose;

% create symbolic variables
[vars,varsDer] = symVariables(sys,true);

% insert symbolic variables into the system equations
[fdyn,fcon,fout] = aux_insertSymVariables(sys,vars);
    
% check if old derivatives can be re-used
[requiredFiles,requiredFiles_out,storedData,deleteAll] = ...
    priv_checkTensorRecomputation(sys,fdyn,fcon,fout,path,options);

% already finished if nothing new required
if ~any([requiredFiles.standard;requiredFiles.int;requiredFiles.higherOrders;...
         requiredFiles_out.standard;requiredFiles_out.int;requiredFiles_out.higherOrders])
    return;
end

% delete all existing files if setttings have changed
if deleteAll
    currDir = pwd;
    if isfolder(path)
        cd(path); delete *.m
    end
    cd(currDir);
end

if options.verbose
    fprintf("Computing multivariate derivatives for dynamic '%s':\n", sys.name)
end

% read value for setting 'simplify' if provided, otherwise choose default
simplifyOpt = getDefaultValue('lagrangeRem.simplify',sys,struct(),options,'options');
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'simplify')
    simplifyOpt = options.lagrangeRem.simplify;
end
paramInt = [];
if isfield(options,'paramInt')
    paramInt = options.paramInt;
end

% generate directory at path
if ~exist(path,'dir')
    mkdir(path);
end
addpath(path);


% --- compute Jacobians ---------------------------------------------------
if options.verbose
    disp('  .. compute symbolic Jacobians');
end

% state/constraint equation
if requiredFiles.standard(1)
    Jdyn = aux_jacobians(fdyn,vars,simplifyOpt);
    Jcon = aux_jacobians(fcon,vars,simplifyOpt);
    [fp,Jp] = aux_parametric(fdyn,Jdyn,vars,paramInt);
    
    % jacobian
    fname = ['jacobian_' sys.name];
    if ~isempty(Jp)
        writeMatrixFile({Jp.x,Jp.u},path,fname,...
            'VarNamesIn',{'x','u','p'},'VarNamesOut',{'A','B'},'BracketSubs',true);
        fname = ['parametricDynamicFile_' sys.name];
        writeMatrixFile({fp},path,fname,...
            'VarNamesIn',{'x','u'},'VarNamesOut',{'f'},'BracketSubs',true);
    elseif ~isempty(vars.y)
        writeMatrixFile({Jdyn.x,Jdyn.u,Jdyn.y,Jcon.x,Jcon.u,Jcon.y},path,fname,...
            'VarNamesIn',{'x','y','u'},'VarNamesOut',{'A','B','C','D','E','F'},'BracketSubs',true);
    else
        writeMatrixFile({Jdyn.x,Jdyn.u},path,fname,...
            'VarNamesIn',{'x','u'},'VarNamesOut',{'A','B'},'BracketSubs',true);
    end
    
    % jacobian_freeParam
    if ~isempty(vars.p)
        fname = ['jacobian_freeParam_' sys.name];
        writeMatrixFile({Jdyn.x,Jdyn.u},path,fname,...
            'VarNamesIn',{'x','u','p'},'VarNamesOut',{'A','B'},'BracketSubs',true);
    end
end
% output equation
if requiredFiles_out.standard(1)
    Jout = aux_jacobians(fout,vars,simplifyOpt);
    [fpout,Jpout] = aux_parametric(fout,Jout,vars,paramInt);

    % out_jacobian
    fname = ['out_jacobian_' sys.name];
    if ~isempty(Jpout)
        writeMatrixFile({Jpout.x,Jpout.u},path,fname,...
            'VarNamesIn',{'x','u','p'},'VarNamesOut',{'C','D'},'BracketSubs',true);
        fname = ['out_parametricDynamicFile_' sys.name];
        writeMatrixFile({fpout},path,fname,...
            'VarNamesIn',{'x','u'},'VarNamesOut',{'f'},'BracketSubs',true);
    elseif ~isempty(vars.y)
        writeMatrixFile({Jout.x,Jout.u},path,fname,...
            'VarNamesIn',{'x','y','u'},'VarNamesOut',{'C','D'},'BracketSubs',true);
    else
        writeMatrixFile({Jout.x,Jout.u},path,fname,...
            'VarNamesIn',{'x','u'},'VarNamesOut',{'C','D'},'BracketSubs',true);
    end
    
    % out_jacobian_freeParam
    if ~isempty(vars.p)
        fname = ['out_jacobian_freeParam_' sys.name];
        writeMatrixFile({Jout.x,Jout.u},path,fname,...
            'VarNamesIn',{'x','u','p'},'VarNamesOut',{'C','D'},'BracketSubs',true);
    end
end


% --- compute Hessians ----------------------------------------------------
if options.verbose
    disp('  .. compute symbolic Hessians');
end
% state equation/constraint equation
if requiredFiles.standard(2) || requiredFiles.int(2)
    J2dyn = aux_hessians(fdyn,vars,simplifyOpt);
    J2con = aux_hessians(fcon,vars,simplifyOpt);
end
% with/without interval arithmetic
if requiredFiles.standard(2)
    fname = ['hessianTensor_' sys.name];
    writeHessianTensorFile(J2dyn,J2con,path,fname,vars,false,options);
end
if requiredFiles.int(2)
    fname = ['hessianTensorInt_' sys.name];
    writeHessianTensorFile(J2dyn,J2con,path,fname,vars,true,options);
end
% output equation
if requiredFiles_out.standard(2) || requiredFiles_out.int(2)
    J2out = aux_hessians(fout,vars,simplifyOpt);
end
% with/without interval arithmetic
if requiredFiles_out.standard(2)
    fname = ['out_hessianTensor_' sys.name];
    writeHessianTensorFile(J2out,[],path,fname,vars,false,options);
end
if requiredFiles_out.int(2)
    fname = ['out_hessianTensorInt_' sys.name];
    writeHessianTensorFile(J2out,[],path,fname,vars,true,options);
end


% --- compute third-order derivatives -------------------------------------
if options.verbose
    disp('  .. compute symbolic third-order derivatives');
end
% state equation/constraint equation
if requiredFiles.standard(3) || requiredFiles.int(3)
    J3dyn = aux_thirdOrderDerivatives(J2dyn,vars,simplifyOpt);
    J3con = aux_thirdOrderDerivatives(J2con,vars,simplifyOpt);
end
% with/without interval arithmetic
if requiredFiles.standard(3)
    fname = ['thirdOrderTensor_' sys.name];
    write3rdOrderTensorFile(J3dyn,J3con,path,fname,vars,false,options);
end
if requiredFiles.int(3)
    fname = ['thirdOrderTensorInt_' sys.name];
    write3rdOrderTensorFile(J3dyn,J3con,path,fname,vars,true,options);
end
% output equation
if requiredFiles_out.standard(3) || requiredFiles_out.int(3)
    J3out = aux_thirdOrderDerivatives(J2out,vars,simplifyOpt);
end
% with/without interval arithmetic
if requiredFiles_out.standard(3)
    fname = ['out_thirdOrderTensor_' sys.name];
    write3rdOrderTensorFile(J3out,[],path,fname,vars,false,options);
end
if requiredFiles_out.int(3)
    fname = ['out_thirdOrderTensorInt_' sys.name];
    write3rdOrderTensorFile(J3out,[],path,fname,vars,true,options);
end


% --- 4th and higher-order tensors (not for output equation) --------------
if any(requiredFiles.higherOrders) % equivalent to options.tensorOrder >= 4
    writeHigherOrderTensorFiles(fdyn,vars,varsDer,path,sys.name,options); 
end


% rehash the folder so that new generated files are used
rehash path;

% save data so that symbolic computations do not have to be re-computed
% d = datetime;
% d.Format = "uuuu-MM-dd_HH-mm_ss";
% filename = sprintf('%s%s%s_%s_%s.mat',path,filesep,sys.name,'_version_',string(d));
filename = sprintf('%s%s%s_%s.mat',path,filesep,sys.name,'lastVersion');
save(filename,'storedData');

% log
if options.verbose
    fprintf("Done.\n")
end

end


% Auxiliary functions -----------------------------------------------------

function options = aux_defaultOptions(sys)
% set required fields
options = struct();

options.tensorOrder = getDefaultValue('tensorOrder',sys,struct(),options,'options');
options.tensorOrderOutput = getDefaultValue('tensorOrderOutput',sys,struct(),options,'options');

options.lagrangeRem.simplify = getDefaultValue('lagrangeRem.simplify',sys,struct(),options,'options');
options.lagrangeRem.tensorParallel = getDefaultValue('lagrangeRem.tensorParallel',sys,struct(),options,'options');
options.lagrangeRem.method = getDefaultValue('lagrangeRem.method',sys,struct(),options,'options');
% no .lagrangeRem.replacements
% no .paramInt

options.verbose = false;

end

function [fdyn,fcon,fout] = aux_insertSymVariables(sys,vars)

if isa(sys,'nonlinearARX')
    fdyn = sys.mFile(vars.x,vars.u);
    fcon = [];
    fout = [];
elseif isempty(vars.y) && isempty(vars.p)
    % class: nonlinearSys, nonlinearSysDT
    fdyn = sys.mFile(vars.x,vars.u);
    fcon = [];
    fout = sys.out_mFile(vars.x,vars.u);
elseif isempty(vars.y) && ~isempty(vars.p)
    % class: nonlinParamSys
    fdyn = sys.mFile(vars.x,vars.u,vars.p);
    fcon = [];
    fout = sys.out_mFile(vars.x,vars.u,vars.p);
elseif ~isempty(vars.y)
    % class: nonlinDASys
    fdyn = sys.dynFile(vars.x,vars.y,vars.u);
    fcon = sys.conFile(vars.x,vars.y,vars.u);
    fout = sys.out_mFile(vars.x,vars.y,vars.u);
end

end

function J = aux_jacobians(f,vars,simplifyOpt)
% jacobians - compute symbolic Jacobians
%
% Inputs:
%    fdyn - symbolic differential equation
%    vars - symbolic variables
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    J - Jacobian

% init
J = struct('x',[],'u',[],'y',[]);

% Jacobians with respect to the state, input, and constraint state
J.x = jacobian(f,vars.x);
J.u = jacobian(f,vars.u);
J.y = jacobian(f,vars.y);

% perform simplification
if strcmp(simplifyOpt,'simplify') 
    J.x = simplify(J.x);
    J.u = simplify(J.u);
    J.y = simplify(J.y);
elseif strcmp(simplifyOpt,'collect')
    J.x = collect(J.x,vars.x);
    J.u = collect(J.u,vars.x);
    J.y = collect(J.y);
end

end

function [fp,Jp] = aux_parametric(f,J,vars,paramInt)
% special function evaluation and derivatives for nonlinear
% dynamic/output equations with parameters linearly influencing the
% derivative

if isempty(vars.p)
    fp = sym([]);
    Jp = struct('x',{},'u',{});
    return
end

% number of states, parameters
numParams = numel(vars.p);

% insert parameters into dynamic file
fp = sym(zeros([size(f,1),1,numParams+1]));
    
% store jacobians with respect to parameters
Jp = struct('x',[],'u',[]);
Jp.x = sym(zeros([size(J.x),numParams+1]));
Jp.u = sym(zeros([size(J.u),numParams+1]));

try
    % part without parameters (-> insert zero for vars.p)
    fp(:,:,1) = subs(f,vars.p,zeros(numParams,1));
    Jp.x(:,:,1) = subs(J.x,vars.p,zeros(numParams,1));
    Jp.u(:,:,1) = subs(J.u,vars.p,zeros(numParams,1));
    % part with parameters
    for i=1:numParams
        fp(:,:,i+1) = subs(f,vars.p,unitvector(i,numParams)) - fp(:,:,1);
        Jp.x(:,:,i+1) = subs(J.x,vars.p,unitvector(i,numParams)) - Jp.x(:,:,1);
        Jp.u(:,:,i+1) = subs(J.u,vars.p,unitvector(i,numParams)) - Jp.u(:,:,1);
    end

    % if parameters are uncertain within an interval
    if isa(paramInt,'interval')
        %normalize
        pCenter = center(paramInt);
        pDelta = rad(paramInt);
    elseif isnumeric(paramInt) && isvector(paramInt)
        pCenter = reshape(paramInt,[],1);
        pDelta = zeros(numel(pCenter),1);
    else
        throw(CORAerror('CORA:notSupported',...
            'paramInt must be an interval or a double vector'));
    end

    for i=1:numParams
        %center
        Jp.x(:,:,1) = Jp.x(:,:,1) + pCenter(i)*Jp.x(:,:,i+1);
        Jp.u(:,:,1) = Jp.u(:,:,1) + pCenter(i)*Jp.u(:,:,i+1);
        %generators
        Jp.x(:,:,i+1) = pDelta(i)*Jp.x(:,:,i+1);
        Jp.u(:,:,i+1) = pDelta(i)*Jp.u(:,:,i+1);
    end

catch
    %%% TODO: find out when this is supposed to occur... and remove?
    Jp = struct('x',{},'u',{});
    disp('Parameters are not linearly influencing the system.');
end

end

function H = aux_hessians(f,vars,simplifyOpt)
% hessians - compute Hessian tensors for differential equation and
%    constraint equation (only nonlinDASys)
%
% Inputs:
%    f - symbolic equation
%    vars - vector of symbolic variables (state, algebraic, input)
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    H - symbolic Hessian tensor

% concatenate symbolic variables (with 'LR')
z = [vars.x;vars.y;vars.u];

% potential simplifications
isSimplify = strcmp(simplifyOpt,'simplify');
isCollect = strcmp(simplifyOpt,'collect');

% Jacobians and Hessians
J_comb = jacobian(f,z);
H = sym('x',[0,0,0]);
for k=1:size(J_comb,1)
    % compute Hessians/2nd-order Jacobians
    H(k,:,:) = jacobian(J_comb(k,:),z);
    if isSimplify
        H(k,:,:) = simplify(H(k,:,:));
    elseif isCollect
        H(k,:,:) = collect(H(k,:,:),z);
    end
end
    
end

function T = aux_thirdOrderDerivatives(H,vars,simplifyOpt)
% thirdOrderDerivatives - compute third-order derivatives for
%    differential equation and constraint equation (only nonlinDASys)
%
% Inputs:
%    H - symbolic Hessian tensor
%    vars - vector of symbolic variables (state, algebraic, input)
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    T - symbolic third-order tensor

[n,nrOfVars,~] = size(H);
T = sym(zeros(n,nrOfVars,nrOfVars,nrOfVars));

% construct vector for which derivative is computed
z = [vars.x;vars.y;vars.u];

% potential simplifications
isSimplify = strcmp(simplifyOpt,'simplify');
isCollect = strcmp(simplifyOpt,'collect');

% compute third-order jacobians using 'LR' variables
for k=1:n
    for l=1:nrOfVars
        % compute 3rd-order Jacobians
        if ~isempty(find(H(k,l,:), 1))
            T(k,l,:,:) = jacobian(reshape(H(k,l,:),[nrOfVars,1]),z);
            if isSimplify
                T(k,l,:,:) = simplify(T(k,l,:,:));
            elseif isCollect
                T(k,l,:,:) = collect(T(k,l,:,:),z);
            end
        end
    end
end
    
end

% ------------------------------ END OF CODE ------------------------------
