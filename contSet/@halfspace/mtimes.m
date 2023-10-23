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
%    M = [0.9 0.2; -0.2 0.9];
%    hs = halfspace([2 1],1);
%    Mhs = M * hs;
% 
%    figure; hold on;
%    plot(hs,[1,2],'b')
%    plot(Mhs,[1,2],'r','FaceAlpha',0.5)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-August-2013
% Last update:   16-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if representsa_(hs,'emptySet',eps)
    return
end

try
    %assume that factor is an invertible matrix
    invMat = inv(M);
    c = invMat.'*hs.c;
    
    if ~all(isfinite(c), 'all')
        % might be due to singular matrix; is checked in catch block
        throw(CORAerror('CORA:wrongValue', 'first', ...
            ['Matrix multiplication lead to an invalid center vector ' ...
            '(center has to be finite).']))
    end

    hs.c = c;
    
catch ME
    
    if diff(size(M)) ~= 0
        throw(CORAerror('CORA:wrongValue','first','square matrix'));
    elseif size(M,2) ~= dim(hs)
        throw(CORAerror('CORA:dimensionMismatch',M,hs));
    elseif abs(det(M)) < eps
        throw(CORAerror('CORA:specialError',...
            'Linear transformation with (near-)singular matrix.'));
    else
        rethrow(ME);
    end
    
end


% ------------------------------ END OF CODE ------------------------------
