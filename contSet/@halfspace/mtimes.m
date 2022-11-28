function hs = mtimes(M,hs)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
%    a halfspace
%
% Syntax:  
%    hs = mtimes(M,hs)
%
% Inputs:
%    M - numerical matrix
%    hs - halfspace object
%
% Outputs:
%    hs - halfspace object
%
% Example: 
%    M = [0.6980 0.7161; -0.7161 0.6980];
%    hs = halfspace([1 1],2);
%    M * hs;
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

% empty case
if isempty(hs)
    return
end

try
    %assume that factor is an invertible matrix
    invMat = inv(M);
    hs.c = invMat.'*hs.c;
    
catch ME
    
    if diff(size(M)) ~= 0
        throw(CORAerror('CORA:wrongValue','first','square matrix'));
    elseif size(M,2) ~= dim(hs)
        throw(CORAerror('CORA:dimensionMismatch',M,hs));
    elseif abs(det(M)) < eps
        throw(CORAerror('CORA:specialError',...
            'Linear transformation with near-singular matrix'));
    else
        rethrow(ME);
    end
    
end


%------------- END OF CODE --------------