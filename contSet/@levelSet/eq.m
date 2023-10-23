function res = eq(ls1,ls2,varargin)
% eq - overloads '==' operator to check if two levelSet objects are equal
%
% Syntax:
%    res = ls1 == ls2
%    res = eq(ls1,ls2)
%    res = eq(ls1,ls2,tol)
%
% Inputs:
%    ls1 - location object
%    ls2 - location object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    % init symbolic variables
%    syms a b c x y z
% 
%    % init level set
%    eq = -x^2 - y^2 + 5;
%    ls1 = levelSet(eq,[x;y],'<=');
% 
%    % different variable names
%    eq = -a^2 - b^2 + 5;
%    ls1_ = levelSet(eq,[a;b],'<=');
% 
%    % compare equal sets
%    ls1 == ls1_
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(ls1,ls2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
