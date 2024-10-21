function nonlinReset = derivatives(nonlinReset,varargin)
% derivatives - compute derivatives of nonlinear reset functions
%
% Syntax:
%    derivatives(nonlinReset)
%    nonlinReset = derivatives(nonlinReset)
%    nonlinReset = derivatives(nonlinReset,fpath)
%    nonlinReset = derivatives(nonlinReset,fpath,fname)
%    nonlinReset = derivatives(nonlinReset,fpath,fname,tensorOrder)
%
% Inputs:
%    nonlinReset - nonlinearReset object
%    fpath - where to store generated files
%    fname - name of reset function, used for file names
%    tensorOrder - tensor order (2 or 3)
%
% Outputs:
%    nonlinReset - nonlinearReset object with set properties for derivatives
%
% Example:
%    nonlinReset = nonlinearReset(@(x,u) [x(1)*x(2); x(2)]);
%    path = [CORAROOT filesep 'models' filesep 'auxiliary'];
%    fname = 'example_derivatives';
%    nonlinReset = derivatives(nonlinReset,path,fname,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       12-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
narginchk(1,5);
defaultPath = [CORAROOT filesep 'models' filesep 'auxiliary' ...
    filesep 'nonlinearReset'];
[fpath,fname,tensorOrder] = ...
    setDefaultValues({defaultPath,'nonlinearReset',2},varargin);

% in derivatives computation
if ~exist(fpath,'dir')
    mkdir(fpath);
end
addpath(fpath);

% note: currently, we override the derivatives at new call of 'derivatives'
% if the file path is the same

% save tensor order
nonlinReset.tensorOrder = tensorOrder;

% generate symbolic variables
x = sym('x',[nonlinReset.preStateDim,1],'real');
u = sym('u',[nonlinReset.inputDim,1],'real');
vars = {x,u};
varNamesIn = {'x','u'};

% Jacobian file
[J_cell,J_han] = derive('FunctionHandle',nonlinReset.f,...
    'Vars',vars,'VarNamesIn',varNamesIn,'VarNamesOut',{'A','B'},...
    'Path',fpath,'FileName',[fname '_jacobian']);
nonlinReset.J = J_han;
% note: when executing J_han(vars{:}), one obtains
%   J_cell{1} = out1, J_cell{2} = out2, ...

if tensorOrder == 1
    return
end

% rewrite for further derivation
J = [J_cell{1}, J_cell{2}];
% J_cell{1}: df1/dx1,dx2,... df2/dx1,dx2,...
% J_cell{2}: df1/du1,du2,... df2/du1,du2,...

% Hessian file
[H_cell,H_han] = derive('SymbolicFunction',J,...
    'Vars',vars,'VarNamesIn',varNamesIn,'VarNamesOut',{'Hx','Hu'},...
    'Path',fpath,'FileName',[fname '_hessian'],...
    'IntervalArithmetic',tensorOrder == 2);
nonlinReset.H = H_han;
% note: when executing H_han(vars{:}), one obtains
%   H_cell{1} = out1, H_cell{2} = out2, ...

if tensorOrder == 2
    return
end

% rewrite for further derivation
H = [H_cell{1}, H_cell{2}];

% third-order tensor file
[~,T_han] = derive('SymbolicFunction',H,...
    'Vars',vars,'VarNamesIn',varNamesIn,'VarNamesOut',{'Tx','Tu'},...
    'Path',fpath,'FileName',[fname '_thirdOrderTensor'],...
    'IntervalArithmetic',tensorOrder == 3);
nonlinReset.T = T_han;

% ------------------------------ END OF CODE ------------------------------
