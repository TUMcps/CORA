function dx =accSysEidFast(x,u)
% accSysEidFast - ???
%
% Syntax:
%    dx =accSysEidFast(x,u)
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

aMax=7;
c1=7.32;

% dx(1,1)=1e-3*x(1)+x(2); %position
% dx(2,1)=1e-3*x(2)+aMax*c1/x(2)*u; %speed

dx(1,1)=x(2); %position
dx(2,1)=aMax*c1/x(2)*u; %speed

% ------------------------------ END OF CODE ------------------------------
