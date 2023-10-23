function res = isequal(hs1,hs2,varargin)
% isequal - checks if two halfspaces are equal
%
% Syntax:
%    res = isequal(hs1,hs2)
%    res = isequal(hs1,hs2,tol)
%
% Inputs:
%    hs1 - halfspace object
%    hs2 - halfspace object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs1 = halfspace([2;4], 4);
%    hs2 = halfspace([-3;5], 3);
%    isequal(hs1,hs2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{hs1,'att','halfspace'};
                {hs2,'att','halfspace'};
                {tol,'att','numeric',{'nonnan','nonnegative','scalar'}}});

% numerical comparison
res = all(withinTol(hs1.c,hs2.c,tol)) && ... % normal vectors
    withinTol(hs1.d,hs2.d,tol); % distances to origin

% ------------------------------ END OF CODE ------------------------------
