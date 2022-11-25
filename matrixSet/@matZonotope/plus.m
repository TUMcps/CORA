function matZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two
% matrix zonotopes or a matrix zonotope with a matrix
%
% Syntax:  
%    [matZ] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope matrix object or numerical matrix
%    summand2 - zonotope matrix object or numerical matrix
%
% Outputs:
%    matZ - matrix zonotpe after Minkowsi addition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%Find a matrix zonotope object
%Is summand1 a matrix zonotope?
if isa(summand1,'matZonotope')
    %initialize resulting zonotope
    matZ=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a matrix zonotope?    
elseif isa(summand2,'matZonotope')
    %initialize resulting zonotope
    matZ=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a zonotope?
if isa(summand,'matZonotope')
    %Calculate minkowski sum
    matZ.center = matZ.center + summand.center;
    
    if isempty(matZ.generator)
        %concatenate matrix generators
        matZ.generator = summand.generator;
    else
        %concatenate matrix generators
        matZ.generator((end+1):(end+summand.gens)) = summand.generator;
    end
    %update number of generators
    matZ.gens=matZ.gens + summand.gens;
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    matZ.center = matZ.center + summand;
end


%------------- END OF CODE --------------