function res = isequal(linsys1,linsys2,varargin)
% isequal - checks if two linear systems are equal
%
% Syntax:
%    res = isequal(linsys1,linsys2)
%    res = isequal(linsys1,linsys2,tol)
%
% Inputs:
%    linsys1 - linearSys object
%    linsys2 - linearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    A = [2 1; -1 2]; B = [1; -1];
%    linsys1 = linearSys(A,B);
%    linsys2 = linearSys(A,B+1e-14);
%    res = isequal(linsys1,linsys2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   30-August-2024 (MW, integrate E and F matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{linsys1,'att','linearSys'};
                {linsys2,'att',{'linearSys','nonlinearSys'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% call nonlinear/isequal if second argument is a nonlinearSys
if isa(linsys2,'nonlinearSys')
    res = isequal(linsys2,linsys1,tol); return
end

% check name
if ~strcmp(linsys1.name,linsys2.name)
    res = false; return
end

% number of states, inputs, outputs
if linsys1.nrOfDims ~= linsys2.nrOfDims || linsys1.nrOfInputs ~= linsys2.nrOfInputs ...
        || linsys1.nrOfOutputs ~= linsys2.nrOfOutputs
    res = false; return
end

% now, sizes of matrices do not need to be checked anymore...

% check system matrix
if ~compareMatrices(linsys1.A,linsys2.A,tol)
    res = false; return
end

% check input matrix
if ~compareMatrices(linsys1.B,linsys2.B,tol)
    res = false; return
end

% check offset (differential equation)
if ~all(withinTol(linsys1.c,linsys2.c,tol))
    res = false; return
end

% check output matrix
if ~all(withinTol(linsys1.C,linsys2.C,tol))
    res = false; return
end

% check feedthrough matrix
if ~all(withinTol(linsys1.D,linsys2.D,tol))
    res = false; return
end

% check offset (output equation)
if ~all(withinTol(linsys1.k,linsys2.k,tol))
    res = false; return
end

% check disturbance (state)
if ~compareMatrices(linsys1.E,linsys2.E,tol)
    res = false; return
end

% check disturbance (output)
if ~compareMatrices(linsys1.F,linsys2.F,tol)
    res = false; return
end

% all checks ok
res = true;

% ------------------------------ END OF CODE ------------------------------
