function r = radius(E,varargin)
% radius - computes the radius of an enclosing hyperball of an ellipsoid
%
% Syntax:
%    r = radius(E)
%    r = radius(E,i)
%
% Inputs:
%    E - ellipsoid object
%    i - number of largest radii to be returned
%
% Outputs:
%    r - radius/vector of radii of enclosing hyperball
%
% Example: 
%    E = ellipsoid([1,0.5;0.5,3],[1;-1]);
%    r = radius(E);
%
%    figure; hold on;
%    plot(E,[1,2],'r');
%    E_circ = ellipsoid(r^2 * eye(2),center(E));
%    plot(E_circ,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Victor Gassmann
% Written:       05-March-2021
% Last update:   19-March-2021 (VG, empty case added)
%                24-March-2022 (VG, change input argument)
%                04-July-2022 (VG, input checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default: only largest radius considered
i = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{E,'att','ellipsoid'};
                {i,'att','numeric',{'integer','positive',@(i) i <= dim(E)}}});

% quick check for empty set
if representsa_(E,'emptySet',eps)
    r = zeros(dim(E),0);
    return
end

% compute eigenvalues
d = eigs(E.Q,i); % since we use Q^{-1} as a shape matrix

% compute radius
r = sqrt(d); % since we use Q^{-1} as a shape matrix

% ------------------------------ END OF CODE ------------------------------
