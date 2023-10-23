function res = isequal(sys1,sys2,varargin)
% isequal - checks if two contDynamics objects are
%
% Syntax:
%    res = isequal(sys1,sys2)
%    res = isequal(sys1,sys2,tol)
%
% Inputs:
%    sys1 - contDynamics object
%    sys2 - contDynamics object
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
% Written:       19-May-2023
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
inputArgsCheck({{sys1,'att','contDynamics'};
                {sys2,'att','contDynamics'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check name
if ~strcmp(sys1.name,sys2.name)
    res = false; return
end

% number of states, inputs, outputs
if sys1.dim ~= sys2.dim || sys1.nrOfInputs ~= sys2.nrOfInputs ...
        || sys1.nrOfOutputs ~= sys2.nrOfOutputs
    res = false; return
end

% all checks ok
res = true;

% ------------------------------ END OF CODE ------------------------------
