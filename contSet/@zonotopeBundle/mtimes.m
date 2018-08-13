function [Zbundle] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
% interval matrix with a zonotope bundle
%
% Syntax:  
%    [Zbundle] = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - zonotope bundle or matrix set or matrix
%    factor2 - zonotope bundle or matrix set or matrix 
%
% Outputs:
%    Zbundle - Zonotope bundle after multiplication
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
%Is factor1 a zonotope bundle?
if isa(factor1,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=factor1;
    %initialize other summand
    factor=factor2;
%Is factor2 a zonotope bundle?    
elseif isa(factor2,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=factor2;
    %initialize other summand
    factor=factor1;  
end


%Calculate multiplication for each zonotope
for i=1:Zbundle.parallelSets
    Zbundle.Z{i}=factor*Zbundle.Z{i};
end


%------------- END OF CODE --------------