function res = isequal(linReset1,linReset2,varargin)
% isequal - checks if two linear reset functions are equal up to some
%    tolerance
%
% Syntax:
%    res = isequal(linReset1,linReset2)
%    res = isequal(linReset1,linReset2,tol)
%
% Inputs:
%    linReset1 - linearReset object
%    linReset2 - linearReset object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    A = eye(2); c1 = [1;-1]; c2 = [1;1];
%    linReset1 = linearReset(A,c1);
%    linReset2 = linearReset(A,c2);
%    isequal(linReset1,linReset1);
%    isequal(linReset1,linReset2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/isequal

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{linReset1,'att','linearReset'};
                {linReset2,'att',{'linearReset','nonlinearReset'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% redirect to nonlinearReset/isequal
if isa(linReset2,'nonlinearReset')
    res = isequal(linReset2,linReset1,tol);
    return
end

% check number of states and inputs in superclass function, this ensures
% that the matrices are of equal size for further checks below
if ~isequal@abstractReset(linReset1,linReset2)
    res = false;
    return
end

% compare matrices A, B, and vector c
res = compareMatrices(linReset1.A,linReset2.A,tol,"equal",true) ...
    && compareMatrices(linReset1.B,linReset2.B,tol,"equal",true) ...
    && compareMatrices(linReset1.c,linReset2.c,tol,"equal",true);

% ------------------------------ END OF CODE ------------------------------
