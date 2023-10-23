function res = isequal(sys1,sys2,varargin)
% isequal - checks if two nonlinear systems are equal
%
% Syntax:
%    res = isequal(sys1,sys2)
%    res = isequal(sys1,sys2,tol)
%
% Inputs:
%    sys1 - nonlinearSys object
%    sys2 - nonlinearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    f = @(x,u) [x(1)^2 - u(1); x(2)];
%    g = @(x,u) [x(2)^2 - u(1); x(1)]; 
%    sys1 = nonlinearSys(f);
%    sys2 = nonlinearSys(g);
%    res = sys1 == sys2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       10-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{sys1,'att','nonlinearSys'};
                {sys2,'att',{'nonlinearSys','linearSys'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% convert linearSys to nonlinearSys
if isa(sys2,'linearSys')
    sys2 = nonlinearSys(sys2);
end

% check name
if ~strcmp(sys1.name,sys2.name)
    res = false; return
end

% quick check: are both output equations are linear/nonlinear
if length(sys1.out_isLinear) ~= length(sys2.out_isLinear) || ...
    any(xor(sys1.out_isLinear,sys2.out_isLinear))
    res = false; return
end

% compare differential equations
if ~aux_compareEquations(sys1.mFile,sys1.dim,sys1.nrOfInputs,sys1.dim,...
        sys2.mFile,sys2.dim,sys2.nrOfInputs,sys2.dim)
    res = false; return
end

% compare output equations
if ~aux_compareEquations(sys1.out_mFile,sys1.dim,sys1.nrOfInputs,sys1.nrOfOutputs,...
        sys2.out_mFile,sys2.dim,sys2.nrOfInputs,sys2.nrOfOutputs)
    res = false; return
end

% all checks ok
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_compareEquations(sys1,x1dim,u1dim,y1dim,sys2,x2dim,u2dim,y2dim)
% sys1, sys2 are function handles which are checked for equality
% sys1: x1dim x u1dim -> y1dim (y1dim = x1dim for differential equations)
% sys2: x2dim x u2dim -> y2dim (y2dim = x2dim for differential equations)

% check number of states, inputs, and outputs
if x1dim ~= x2dim || y1dim ~= y2dim
    res = false; return
end

% inputs can also be dummy inputs, so length does not necessarily have to
% be the same for the equations to be equal -> choose max
udim = max([u1dim,u2dim]);

% instantiate symbolic variables to evaluate nonlinear equations
x_sym = sym('x',[x1dim,1],'real');
u_sym = sym('u',[udim,1],'real');

% evaluate equations
f1 = sys1(x_sym,u_sym);
f2 = sys2(x_sym,u_sym);

% 1. symbolic difference should be zero
sym_diff = simplify(f1-f2);
idx_diff_zero = logical(sym_diff == 0);
if all(idx_diff_zero)
    res = true; return
end


% 2. plug in ('sufficiently many') values

% number of test evaluations
nrTests = 1000;
for i=1:nrTests

    % random values for state and input
    x_rand = randn(x1dim,1);
    u_rand = randn(udim,1);

    % plug in both function handles
    f1_rand = sys1(x_rand,u_rand);
    f2_rand = sys2(x_rand,u_rand);

    % numerical comparison
    comp = withinTol(f1_rand,f2_rand);
    % only check dims whose difference was not symbolically already zero
    if ~all(comp(~idx_diff_zero))
        % counter-example found
        res = false; return
    end
end

end

% ------------------------------ END OF CODE ------------------------------
