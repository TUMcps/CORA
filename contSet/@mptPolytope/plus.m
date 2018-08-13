function [P] = plus(summand1,summand2)
% plus - overloaded '+' operator for the addition of a vector to a
% mptPolytope; Minkowski addition is also supported
%
% Syntax:  
%    [P] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - mptPolytope object or numerical vector
%    summand2 - mptPolytope object or numerical vector
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
% Written:      20-October-2010
% Last update:  24-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Find a pplPolytope object
%Is summand1 a mptPolytope?
if isa(summand1,'mptPolytope')
    %initialize resulting mptPolytope
    P=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a mptPolytope?    
elseif isa(summand2,'mptPolytope')
    %initialize resulting mptPolytope
    P=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a vector?
if isnumeric(summand)
    try %MPT V3
        P.P = P.P + summand;
    catch %MPT V2
        dim = length(summand);
        P.P = range(P.P, eye(dim), summand);
    end
    
%is summand a polytope?
elseif isa(summand,'mptPolytope')
    P.P = P.P + summand.P;
    
%something else?    
else
    P=[];
    disp('this operation is not implemented');
end

%------------- END OF CODE --------------