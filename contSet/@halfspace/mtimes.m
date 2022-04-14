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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-August-2013
% Last update:  16-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

try
    %assume that factor is an invertible matrix
    invMat = inv(factor);
    h.c = invMat.'*h.c;
    
catch ME
    
    if isempty(h)
        % empty halfspace
        [msg,id] = errEmptySet();
        error(id,msg);
    elseif abs(det(factor)) < eps
        error("Linear transformation with near-singular matrix");
    else
        rethrow(ME);
    end
    
end


%------------- END OF CODE --------------