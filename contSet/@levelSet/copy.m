function ls_out = copy(ls)
% copy - copies the levelSet object (used for dynamic dispatch)
%
% Syntax:
%    ls_out = copy(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    ls_out - copied levelSet object
%
% Example: 
%    syms x y;
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    ls_out = copy(ls);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
ls_out = levelSet(ls);

% ------------------------------ END OF CODE ------------------------------
