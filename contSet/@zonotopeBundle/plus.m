function [Zbundle] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a zonotope
% bundle with a zonotope or with a vector
%
% Syntax:  
%    [Zbundle] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope bundle or zonotope or numerical vector
%    summand2 - zonotope bundle or zonotope or numerical vector
%
% Outputs:
%    Zbundle - Zonotpe bundle after Minkowsi addition
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
% Written:      09-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope bundle object
%Is summand1 a zonotope bundle?
if isa(summand1,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a zonotope bundle?    
elseif isa(summand2,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=summand2;
    %initialize other summand
    summand=summand1;  
end

% over-approximate the set if the summand is a zontope bundle
if isa(summand,'zonotopeBundle')
   summand = zonotope(interval(summand)); 
end

%Calculate minkowski sum for each zonotope
for i=1:Zbundle.parallelSets
    Zbundle.Z{i}=Zbundle.Z{i}+summand;
end


%------------- END OF CODE --------------