function n = dim(E)
% dim - returns the dimension of the ambient space of an ellipsoid
%
% Syntax:
%    n = dim(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    n = dim(E) 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       15-September-2019 
% Last update:   16-March-2021 (comp independent of property)
%                04-July-2022 (VG, support class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = zeros(size(E));
for i=1:numel(E)
    try
        n(i) = length(E(i).q);
    catch ME
        if isemptyobject(E(i))
            n(i) = 0;
        else
            rethrow(ME);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
