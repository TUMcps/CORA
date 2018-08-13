function [V] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
% interval matrix with vertices
%
% Syntax:  
%    [V] = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - matrix or vertices object 
%    factor2 - matrix or vertices object 
%
% Outputs:
%    V - Vertices after multiplication with a matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author: Matthias Althoff
% Written: 07-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%Find a vertices object
%Is factor1 a vertices object?
if strcmp('vertices',class(factor1))
    %initialize resulting zonotope
    V=factor1;
    %initialize other summand
    matrix=factor2;
    %compute result
    V.V=V.V*matrix;  
%Is factor2 a vertices object?    
elseif strcmp('vertices',class(factor2))
    %initialize resulting zonotope
    V=factor2;
    %initialize other summand
    matrix=factor1;  
    %compute result
    V.V=matrix*V.V;     
end


%------------- END OF CODE --------------