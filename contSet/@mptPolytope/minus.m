function [P] = minus(minuend,subtrahend)
% minus - overloaded '-' operator for the subtraction of a vector from an
% mptPolytope; Minkowski difference is also supported
%
% Syntax:  
%    [P] = minus(summand1,summand2)
%
% Inputs:
%    minuend - mptPolytope object or numerical vector
%    subtrahend - mptPolytope object or numerical vector
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      22-July-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%Is subtrahend a mptPolytope?
if isa(subtrahend,'mptPolytope')
    P = minuend;
    P.P = minuend.P - subtrahend.P;

%is subtrahend a vector?
elseif isnumeric(subtrahend)
    P = minuend + (-subtrahend);
end
    

%------------- END OF CODE --------------