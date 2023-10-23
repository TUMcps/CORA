function dirs = randEqdistDirections(dim,sections)
% randEqdistDirections - computes evenly distributed directions in all
%    dimensions
%
% Syntax:
%    dirs = randEqdistDirections(dim,sections)
%
% Inputs:
%    dim - dimension
%    sections - number of sections in each dimension
%
% Outputs:
%    dirs - matrix of directions
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-September-2012
% Last update:   04-May-2020 (MW, rename from "additionalGenerators")
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%first unit vector
dirs = [1;zeros(dim-1,1)];

%increment
incr = pi/(sections+1);

%loop over all dimensions
for iDim = 1 : dim-1
    %set up Gnew
    G = dirs;
    dirs = [];
    for iAngle = 1 : sections
        %determine angle
        angle = iAngle*incr;
        
        %construct rotational matrix
        RotMat = eye(dim);
        RotMat(iDim, iDim) = cos(angle);
        RotMat(iDim, iDim+1) = -sin(angle);
        RotMat(iDim+1, iDim) = sin(angle);
        RotMat(iDim+1, iDim+1) = cos(angle);
        
        %compute new generators
        for iGen = 1 : length(G(1,:))
            dirs(:,end+1) = RotMat*G(:,iGen);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
