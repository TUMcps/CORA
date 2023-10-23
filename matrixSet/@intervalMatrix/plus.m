function intMat = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two interval
%    matrices or a interval matrix with a matrix
%
% Syntax:
%    intMat = plus(summand1,summand2)
%
% Inputs:
%    summand1 - interval matrix object or numerical matrix
%    summand2 - interval matrix object or numerical matrix
%
% Outputs:
%    intMat - interval matrix after Minkowsi addition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   05-August-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find a matrix zonotope object
%Is summand1 a matrix zonotope?
if isa(summand1,'intervalMatrix')
    %initialize resulting zonotope
    intMat=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a matrix zonotope?    
elseif isa(summand2,'intervalMatrix')
    %initialize resulting zonotope
    intMat=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand an interval matrix?
if isa(summand,'intervalMatrix')
    %Calculate minkowski sum
    intMat.int = intMat.int + summand.int;

%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    intMat.int = intMat.int + summand;
end

% ------------------------------ END OF CODE ------------------------------
