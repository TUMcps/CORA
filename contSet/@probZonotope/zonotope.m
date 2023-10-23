function Z = zonotope(probZ,m)
% zonotope - converts a probabilistic zonotope to a common zonotope where for
%    each generator, a m-sigma interval is taken
%
% Syntax:
%    Z = zonotope(probZ,m)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    m - integer determining the size of the intervals
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    Z = zonotope(probZ);
% 
%    figure; hold on;
%    plot(probZ); plot(Z,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       03-August-2007 
% Last update:   26-February-2008
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin==1
    m=probZ.gamma;
end

%reduce probabilistic zonotope
pZred=probReduce(probZ);

%new generators
newG=m*pZred.g;

%build new zonotope
c = center(probZ);
G = [probZ.Z(:,2:end),newG];
Z = zonotope(c,G);

% ------------------------------ END OF CODE ------------------------------
