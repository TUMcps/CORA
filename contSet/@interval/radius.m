function r = radius(I)
% radius - Computes radius of enclosing hyperball of an interval 
%
% Syntax:
%    r = radius(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    r - radius of enclosing hyperball
%
% Example: 
%    I = interval([-2;1],[4;3]);
%    r = radius(I);
%
%    figure; hold on;
%    plot(I,[1,2],'r');
%    E = ellipsoid(r^2 * eye(2),center(I));
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       22-July-2016
% Last update:   02-September-2019 (rename enclosingRadius -> radius)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

r = sqrt(sum(rad(I).^2));

% ------------------------------ END OF CODE ------------------------------
