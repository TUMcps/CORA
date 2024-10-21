function res = isequal(nonlinReset1,nonlinReset2,varargin)
% isequal - checks if two nonlinear reset functions, or a nonlinear reset
%    function and a linear reset function, are equal up to some tolerance
%
% Syntax:
%    res = isequal(nonlinReset1,nonlinReset2)
%    res = isequal(nonlinReset1,nonlinReset2,tol)
%
% Inputs:
%    nonlinReset1 - nonlinearReset object, linearReset object
%    nonlinReset2 - nonlinearReset object, linearReset object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    f = @(x,u) [-x(1)*x(2); x(2) - u(1)];
%    g = @(x,u) [-x(1)*x(2); x(1) - u(1)];
%    nonlinReset1 = nonlinearReset(f);
%    nonlinReset2 = nonlinearReset(g);
%    isequal(nonlinReset1,nonlinReset1);
%    isequal(nonlinReset1,nonlinReset2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearReset/isequal, isequalFunctionHandle

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{nonlinReset1,'att','nonlinearReset'};
                {nonlinReset2,'att',{'nonlinearReset','linearReset'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check number of states and inputs in superclass function
if ~isequal@abstractReset(nonlinReset1,nonlinReset2)
    res = false;
    return
end

% convert to second input argument to nonlinearReset object
if isa(nonlinReset2,'linearReset')
    nonlinReset2 = nonlinearReset(nonlinReset2);
end

% compare function handles
res = isequalFunctionHandle(nonlinReset1.f,nonlinReset2.f);

% ------------------------------ END OF CODE ------------------------------
