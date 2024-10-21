function nonlinReset_ = lift(nonlinReset,N,M,stateBind,inputBind,id)
% lift - lifts a nonlinear reset function to a higher-dimensional space
%
% Syntax:
%    nonlinReset_ = lift(linReset,N,M,stateBind,inputBind,id)
%
% Inputs:
%    linReset - nonlinearReset object
%    N - dimension of the higher-dimensional state space
%    M - dimension of higher-dimensional input space
%    stateBind - states of the high-dimensional space that correspond to
%                the states of the low-dimensional reset object
%    inputBind - inputs of the high-dimensional space that correspond to
%                the inputs of the low-dimensional reset object
%    id - true/false whether identity reset function should be used for all
%         other states
%
% Outputs:
%    nonlinReset_ - lifted nonlinearReset object
%
% Example:
%    f = @(x,u) [x(1)*x(2); x(2)];
%    nonlinReset = nonlinearReset(f);
%    N = 6; M = 5; stateBind = [2,3]; inputBind = 3; id = false;
%    nonlinReset_ = lift(nonlinReset,N,M,stateBind,inputBind,id);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearReset/lift

% Authors:       Mark Wetzlinger
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% projection only allowed for R^n -> R^n (and no inputs)
if nonlinReset.preStateDim ~= nonlinReset.postStateDim
    throw(CORAerror('CORA:notSupported',...
        ['Projection of reset functions to higher-dimensional spaces '...
        'only supported for R^n -> R^n.']));
end

% read out old number of states and inputs
n_old = nonlinReset.preStateDim;
m_old = nonlinReset.inputDim;

% init symbolic variables
x_old = sym('x',[n_old,1],'real');
u_old = sym('u',[m_old,1],'real');

% new symbolic variables of larger length
% (note: we require size+1 for correct indexing via matlabFunction below)
x_new = sym('x',[N+1,1],'real');
u_new = sym('u',[M+1,1],'real');

% init new function of correct size
if id
    f_sym = x_new(1:end-1);
else
    f_sym = sym(zeros(N,1));
end

% evaluate nonlinear function with old symbolic variables and replace by
% new ones in lifted space
f_sym_old = nonlinReset.f(x_old,u_old);
f_sym(stateBind) = ...
    subs(f_sym_old,[x_old;u_old],[x_new(stateBind);u_new(inputBind)]);

% init function handle from symbolic
fproj = matlabFunction(f_sym,'Vars',{x_new,u_new});

% set prestate/input dimensions to correct value (they are determined
% automatically by the constructor, but will only be the maximum occuring
% indicies in f_sym, not necessarily N or M)
nonlinReset_ = nonlinearReset(fproj,N,M,size(f_sym,1));

% ------------------------------ END OF CODE ------------------------------
