function G = generators(E)
% generators - returns the generator matrix of an ellipsoid in generator
%    representation
%
% Syntax:
%    G = generators(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    G - generator matrix
%
% Example: 
%    E = ellipsoid([1 1 0; 1 2 0; 0 0 0]);
%    G = generators(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       05-March-2024 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Q = E.Q;
if isempty(Q)
    G = zeros([dim(E) 0]);
    return
end


[U,D,~] = svd(Q);

G = U * sqrt(D);

% ------------------------------ END OF CODE ------------------------------
