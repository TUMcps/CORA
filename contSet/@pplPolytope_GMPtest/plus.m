function [P] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the addition of a vector to a
% pplPolytope; Minkowski addition is not supported
%
% Syntax:  
%    [P] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - pplPolytope object or numerical vector
%    summand2 - pplPolytope object or numerical vector
%
% Outputs:
%    P - pplPolytope object
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
% Written:      20-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%Find a pplPolytope object
%Is summand1 a pplPolytope?
if strcmp('pplPolytope',class(summand1))
    %initialize resulting pplPolytope
    P=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a pplPolytope?    
elseif strcmp('pplPolytope',class(summand2))
    %initialize resulting pplPolytope
    Z=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a vector?
if isnumeric(summand)
    P.d=P.d+P.C*summand;
end

%------------- END OF CODE --------------