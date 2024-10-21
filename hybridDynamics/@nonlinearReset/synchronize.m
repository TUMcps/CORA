function nonlinReset_sync = synchronize(nonlinResets,varargin)
% synchronize - synchronize nonlinear reset functions of equal pre-state
%    and post-state dimensions by concatenating their respective
%    nonlinear functions, which corresponds to the reset functions being
%    evaluated synchronously
%
% Syntax:
%    nonlinReset_sync = synchronize(nonlinResets)
%    nonlinReset_sync = synchronize(nonlinResets,idStates)
%
% Inputs:
%    nonlinReset - class array of nonlinearReset objects
%    idStates - indices of states that are mapped by identity
%
% Outputs:
%    nonlinReset_sync - nonlinearReset object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearReset/synchronize

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       04-April-2022
% Last update:   01-July-2022
%                14-January-2023 (MW, handle states unaffected by sync)
% Last revision: 13-October-2024 (MW, moved from transition class)

% ------------------------------ BEGIN CODE -------------------------------

% set default value for identity mapping
preStateDim = nonlinResets(1).preStateDim;
postStateDim = nonlinResets(1).postStateDim;
idStates = setDefaultValues({[]},varargin);

% all reset functions need to have the same mapping
preStateDims = arrayfun(@(r) r.preStateDim,nonlinResets,'UniformOutput',true);
postStateDims = arrayfun(@(r) r.postStateDim,nonlinResets,'UniformOutput',true);
if any(postStateDims ~= postStateDims(1)) || any(preStateDims ~= preStateDims(1))
    throw(CORAerror('CORA:wrongValue','first',...
        'All reset functions must map the same spaces: R^n -> R^m.'));
end
inputDim = nonlinResets(1).inputDim;

% synchronize nonlinear resets (again: hack for matlabFunction below)
x_sym = sym('x',[preStateDim+1,1],'real');
u_sym = sym('u',[inputDim+1,1],'real');

% init symbolic function using identity mappings (note that we use the same
% symbolic variable 'x' here but with the size of the output vector)
x_sym_out = sym('x',[postStateDim,1],'real');
f_sym = sym(zeros(postStateDim,1));
f_sym(idStates) = x_sym_out(idStates);

% loop over all nonlinear reset functions, add them up
for i=1:numel(nonlinResets)
    f_sym = f_sym + nonlinResets(i).f(x_sym(1:end-1),u_sym(1:end-1));
end
% convert to function handle
f_sync = matlabFunction(f_sym,'Vars',{x_sym,u_sym});

% initialize synchronized nonlinear reset function (ensure that correct
% dimensions are set because highest indices do not necessarily occur)
nonlinReset_sync = nonlinearReset(f_sync,preStateDim,inputDim,size(f_sym,1));

% ------------------------------ END OF CODE ------------------------------
