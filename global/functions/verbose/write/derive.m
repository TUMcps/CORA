function [symDerivative,handle] = derive(varargin)
% derive - compute customizable derivatives of nonlinear functions and save
%    them as m-files; the number of input arguments is the length of the
%    cell array for the provided set of variables (this can also be read
%    automatically from the provided function handle), the number of output
%    arguments is the derived tensor with respect to each variable in the
%    set of variables
%
% Syntax:
%    derive('FunctionHandle',f)
%    derive('SymbolicFunction',f_sym)
%    ...
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'FunctionHandle',f> - function handle
%       <'SymbolicFunction',f_sym> - symbolic function (nD array of sym)
%       <'Vars',vars> - set of variables w.r.t which function is derived
%       <'VarNamesIn',varNamesIn> - variable names for input arguments
%       <'VarNamesOut',varNamesOut> - variable names for output arguments
%       <'Path',fpath> - path where to save generated m-file
%       <'FileName',fname> - name of generated file
%       <'Verbose',verbose> - true/false for verbose output
%       <'IntervalArithmetic',isInt> - true/false whether derivative is to
%           be evaluated using interval arithmetic
%       (extend to .lagrangeRem, parametric, ...)
%
% Outputs:
%    symDerivative - sym object containing the symbolic derivative
%    handle - handle to the generated m-file
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: writeMatrixFile

% Authors:       Mark Wetzlinger
% Written:       12-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
end

% read input arguments
NVpairs = varargin(1:end);
% check list of name-value pairs
checkNameValuePairs(NVpairs,{'FunctionHandle','SymbolicFunction',...
    'Vars','VarNamesIn','VarNamesOut','Path','FileName','IntervalArithmetic','Verbose'});
% function handle given?
[NVpairs,f] = readNameValuePair(NVpairs,'FunctionHandle',...
    @(x) isa(x,'function_handle'),@()[]);
% symbolic function given?
[NVpairs,f_sym] = readNameValuePair(NVpairs,'SymbolicFunction',...
    @(x) isa(x,'sym'),sym('f',[0,0]));
% path for storage given?
defaultPath = [CORAROOT filesep 'models' filesep 'auxiliary'];
[NVpairs,fpath] = readNameValuePair(NVpairs,'Path',@ischar,defaultPath);
% name of system given?
[NVpairs,fname] = readNameValuePair(NVpairs,'FileName',@ischar,'derivative');
% interval arithmetic true/false?
[NVpairs,intervalOn] = readNameValuePair(NVpairs,'IntervalArithmetic',@islogical,false);
% verbose output?
[NVpairs,verbose] = readNameValuePair(NVpairs,'Verbose',@islogical,false);

% symbolic variables, names for input/output arguments given?
% ...their default values are a bit more complicated
[NVpairs,vars] = readNameValuePair(NVpairs,'Vars');
[NVpairs,varNamesIn] = readNameValuePair(NVpairs,'VarNamesIn');
[NVpairs,varNamesOut] = readNameValuePair(NVpairs,'VarNamesOut');

% set all remaining values
[vars,varNamesIn,varNamesOut] = aux_setDefaultValues(vars,varNamesIn,varNamesOut);

% check values
inputArgsCheck({{f,'att','function_handle','scalar'};...
                {f_sym,'att','sym'};...
                {vars,'att','cell'};...
                {varNamesIn,'att','cell'};...
                {varNamesOut,'att','cell'};...
                {fpath,'att','char'};...
                {fname,'att','char'};...
                {intervalOn,'att','logical','scalar'};...
                {verbose,'att','logical','scalar'}});

% evaluate function handle if necessary, otherwise use symbolic function
f_sym = aux_getSymbolicFunction(f_sym,f,vars);

% derive using built-in 'jacobian' function (from symbolic toolbox)
symDerivative = aux_derive(f_sym,vars);

% skip file generation if path is 'none'
if strcmp(fpath,'none')
    handle = [];
    return
end

% substitute x1 -> xL1R, etc. so that we can write the file correctly
numVars = numel(vars);
vars_LR = arrayfun(@(i) sym([varNamesIn{i} 'L%dR'],[numel(vars{i})+1,1],'real'),...
    1:numVars,'UniformOutput',false);
vars_LR = cellfun(@(sym_array) sym_array(1:end-1),vars_LR,'UniformOutput',false);
oldVars = vertcat(vars{:});
newVars = vertcat(vars_LR{:});
symDerivative_LR = cellfun(@(M) subs(M,oldVars,newVars),...
    symDerivative,'UniformOutput',false);

% generate file
handle = writeMatrixFile(symDerivative_LR,fpath,fname,...
    'varNamesIn',varNamesIn,'varNamesOut',varNamesOut,...
    'BracketSubs',true,'IntervalArithmetic',intervalOn);

end


% Auxiliary functions -----------------------------------------------------

function symDerivative = aux_derive(f_sym,vars)

% check size of symbolic expression to call 'jacobian' correctly
sz = size(f_sym);
% number of variables
numVars = numel(vars);
symDerivative = cell(numVars,1);

% support only up until 3D for now... nD version probably requires
% recursion (no time one week before v2025)
if numel(sz) > 3
    throw(CORAerror('CORA:notSupported',...
        'The function derive currently only supports max. 3D sym matrices.'));
end

if isvector(f_sym)
    % vectors are ok as they are... reshape has no effect on 'jacobian' but
    % might make it more obvious what is happening
    f_sym = reshape(f_sym,[],1);
    symDerivative = cellfun(@(var) jacobian(f_sym,var),vars,...
        'UniformOutput',false);

elseif numel(sz) == 2
    % 2D array (no vector) -> 3D tensor
    for i=1:numVars
        for j=1:sz(1)
            symDerivative{i}(j,:,:) = jacobian(f_sym(j,:),vars{i});
        end
    end

elseif numel(sz) == 3
    % 3D array -> 4D tensor
    for i=1:numVars
        for k=1:sz(1)
            for l=1:sz(2)
                symDerivative{i}(k,l,:,:) = ...
                    jacobian(reshape(f_sym(k,l,:),[sz(2),1]),vars{i});
            end
        end
    end

end

end

function f_sym = aux_getSymbolicFunction(f_sym,f,vars)

% prefer symbolic function over function handle (warning if both given)
% note: we require the string because the default value (f = @() [];) does
% not return true when inserted into isempty)
f_text = replace(func2str(f),{'@','(',')','[',']'},'');
if ~isempty(f_text) && ~isempty(f_sym)
    CORAwarning("CORA:global",...
        "Either provide function handle or symbolic function.\nSymbolic function used.");
end

% evaluate symbolic function
if isempty(f_sym)
    try
        f_sym = f(vars{:});
    catch ME
        throw(CORAerror('CORA:specialError',...
            'Variables do not match given function handle.'));
    end
end

end

function [vars,varNamesIn,varNamesOut] = ...
    aux_setDefaultValues(vars,varNamesIn,varNamesOut)

% if vars not given, use inputArgsLength to generate variables of
% appriopriate length: example f(x,u) with x in R^3, u in R^2 yields
%    vars{1} = [in1_1; in1_2];
%    vars{2} = [in2_1; in2_2; in2_3];
% (we use 'in' analogously to matlabFunction)
if isempty(vars)
    count = inputArgsLength(f);
    numInputArgs = numel(count);
    maxInputArgLength = max(count);
    vars_raw = mat2cell(sym('in',[numInputArgs,maxInputArgLength]).',...
        maxInputArgLength,ones(1,numInputArgs));
    vars = arrayfun(@(ii) vars_raw{ii}(1:count(ii)),1:numInputArgs,...
        'UniformOutput',false);
end

% default input arguments
if isempty(varNamesIn)
    varNamesIn = arrayfun(@(i) sprintf('in%i',i),1:numel(vars),...
        'UniformOutput',false);
end

% default: variable names out1, out2, ...
if isempty(varNamesOut)
    varNamesOut = arrayfun(@(i) sprintf('out%i',i),1:numel(vars),...
        'UniformOutput',false);
end

end

% ------------------------------ END OF CODE ------------------------------
