function [h] = mtimes(factor,h)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
% a halfspace
%
% Syntax:  
%    [h] = mtimes(factor,h)
%
% Inputs:
%    factor - numerical matrix
%    h - halfspace object
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
% See also: plus

% Author:       Matthias Althoff
% Written:      26-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%assume that factor is an invertible matrix
invMat = inv(factor);
h.c = invMat.'*h.c;


%------------- END OF CODE --------------