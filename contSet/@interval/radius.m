function r = radius(obj)
% radius - Computes radius of enclosing hyperball of an interval 
%
% Syntax:  
%    r = radius(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    r - radius of enclosing hyperball
%
% Example: 
%    I = interval.generateRandom(2);
%    r = radius(I);
%
%    figure
%    hold on
%    plot(I,[1,2],'r');
%    E = ellipsoid(r^2 * eye(2),center(I));
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      22-July-2016
% Last update:  02-Sep-2019 (rename enclosingRadius -> radius)
% Last revision:---

%------------- BEGIN CODE --------------

%compute radius
r = sqrt(sum(rad(obj).^2));

%------------- END OF CODE --------------