function res = ne(linReset1,linReset2,varargin)
% ne - checks if two linear reset functions are equal up to some tolerance
%
% Syntax:
%    res = linReset ~= linReset2
%    res = ne(linReset1,linReset2)
%    res = ne(linReset1,linReset2,tol)
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
%    linReset1 ~= linReset1;
%    linReset1 ~= linReset2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearReset/isequal

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% redirect to isequal
res = ~isequal(linReset1,linReset2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
