function res = isequal(sys1,sys2,varargin)
% isequal - checks if two linear systems are equal
%
% Syntax:
%    res = isequal(sys1,sys2)
%    res = isequal(sys1,sys2,tol)
%
% Inputs:
%    sys1 - linearSys object
%    sys2 - linearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    A = [2 1; -1 2]; B = [1; -1];
%    sys1 = linearSys(A,B);
%    sys2 = linearSys(A,B+1e-14);
%    res = isequal(sys1,sys2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
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
inputArgsCheck({{sys1,'att','linearSys'};
                {sys2,'att',{'linearSys','nonlinearSys'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% call nonlinear/isequal if second argument is a nonlinearSys
if isa(sys2,'nonlinearSys')
    res = isequal(sys2,sys1,tol); return
end

% check name
if ~strcmp(sys1.name,sys2.name)
    res = false; return
end

% number of states, inputs, outputs
if sys1.dim ~= sys2.dim || sys1.nrOfInputs ~= sys2.nrOfInputs ...
        || sys1.nrOfOutputs ~= sys2.nrOfOutputs
    res = false; return
end

% now, sizes of matrices do not need to be checked anymore...

% check system matrix
if ~all(all(withinTol(sys1.A,sys2.A,tol)))
    res = false; return
end

% check input matrix
if ~all(all(withinTol(sys1.B,sys2.B,tol)))
    res = false; return
end

% check offset (differential equation): requires additional check for
% emptiness as sys.c can be empty...
if xor(isempty(sys1.c),isempty(sys2.c))
    res = false; return
elseif ~all(withinTol(sys1.c,sys2.c,tol))
    res = false; return
end

% check output matrix
if ~all(withinTol(sys1.C,sys2.C,tol))
    res = false; return
end

% check feedthrough matrix
if ~all(withinTol(sys1.D,sys2.D,tol))
    res = false; return
end

% check offset (output equation)
if ~all(withinTol(sys1.k,sys2.k,tol))
    res = false; return
end

% all checks ok
res = true;

% ------------------------------ END OF CODE ------------------------------
