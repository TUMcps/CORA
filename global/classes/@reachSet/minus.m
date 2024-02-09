function R = minus(R1,R2)
% minus - Overloaded '-' operator for the Minkowski addition of a
%     set/vector and the negation of a second set/vector
%
% Syntax:
%    R = minus(R,S)
%
% Inputs:
%    R1 - numeric, contSet, or reachSet object 
%    R2 - numeric, contSet, or reachSet object 
%
% Outputs:
%    R - resulting tranformed reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/minus

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

R = R1 + (-R2);

% ------------------------------ END OF CODE ------------------------------
