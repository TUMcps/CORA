function [dx]=accSysEidFastBicycle(x,u)
% accSysEidFastBicycle - ???
%
% Syntax:
%    [dx]=accSysEidFastBicycle(x,u)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       ???
% Written:       ???
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

aMax=2.5;
c1=2.5;

% dx(1,1)=1e-3*x(1)+x(2); %position
% dx(2,1)=1e-3*x(2)+aMax*c1/x(2)*u; %speed

dx(1,1)=x(2); %position
dx(2,1)=aMax*c1/x(2)*u; %speed

% ------------------------------ END OF CODE ------------------------------
