function C = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with 
% a capsule
%
% Syntax:  
%    C = mtimes(matrix,C)
%
% Inputs:
%    matrix - numerical matrix
%    C - capsule object 
%
% Outputs:
%    C - capsule after multiplication with a matrix
%
% Example: 
%    C = capsule([1; 1], [0; 1], 0.5);
%    matrix=[0 1; 1 0];
%    plot(C);
%    hold on
%    C = matrix*C;
%    plot(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

%Find a capsule object
%Is factor1 a capsule?
if isa(factor1,'capsule')
    %initialize resulting zonotope
    C = factor1;
    %initialize other summand
    matrix = factor2;
%Is factor2 a zonotope?    
elseif isa(factor2,'capsule')
    %initialize resulting zonotope
    C = factor2;
    %initialize other summand
    matrix = factor1;  
end

%numeric matrix
if isnumeric(matrix)
    % new center
    C.c = matrix*C.c;
    % new generator
    C.g = matrix*C.g;
    % new axes of ellpsoid of transformed ball
    newAxes = eig(matrix*matrix');
    C.r = C.r*max(newAxes);
    
%something else?    
else
    % throw error for given arguments
    error(noops(factor1,factor2));
end

%------------- END OF CODE --------------