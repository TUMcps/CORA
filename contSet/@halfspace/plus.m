function [h] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the addition of a vector with a
% halfspace
%
% Syntax:  
%    [h] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - halfspace object or numerical vector
%    summand2 - halfspace object or numerical vector
%
% Outputs:
%    h - halfspace object
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
% Written:      28-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%Find a halfspace object
%Is summand1 a halfspace?
if isa(summand1,'halfspace')
    %initialize resulting zonotope
    h=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a zonotope?    
elseif isa(summand2,'halfspace')
    %initialize resulting zonotope
    h=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a zonotope?
if isnumeric(summand)
    %Calculate minkowski sum
    h.d = h.d + h.c.'*summand;
else
    h = [];
end

%------------- END OF CODE --------------