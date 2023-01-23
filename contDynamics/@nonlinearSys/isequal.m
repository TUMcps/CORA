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

% Author:       Mark Wetzlinger, Victor Gassmann
% Written:      10-January-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
if ~aux_compareEquations(sys1.mFile,sys2.mFile)
    res = false; return
end

% compare output equations
if ~aux_compareEquations(sys1.out_mFile,sys2.out_mFile)
    res = false; return
end

% all checks ok
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_compareEquations(sys1,sys2)
% sys1, sys2 are function handles which are checked for equality
% xu(1) is the number of states, xu(2) is the number of inputs

% obtain number of states, inputs, and outputs
[xu,y] = inputArgsLength(sys1);
[xu_,y_] = inputArgsLength(sys2);

% check number of states, inputs, and outputs
if any(xu ~= xu_) || y ~= y_
    res = false; return
end

% instantiate symbolic variables to evaluate nonlinear equations
x_sym = sym('x',[xu(1),1],'real');
u_sym = sym('u',[xu(2),1],'real');

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
    x_rand = randn(xu(1),1);
    u_rand = randn(xu(2),1);

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

%------------- END OF CODE --------------