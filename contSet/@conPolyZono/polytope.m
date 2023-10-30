function P = polytope(cPZ,varargin)
% polytope - Encloses a constrained polynomial zonotope by a polytope
%
% Syntax:
%    P = polytope(cPZ)
%    P = polytope(cPZ,type)
%    P = polytope(cPZ,type,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    type - algorithm used to compute polytope enclosure ('linearize',
%           'extend', or 'all')
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    P - polytope object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    EC = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    E = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    P = polytope(cPZ);
%
%    figure; hold on
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(P,[1,2],'b');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope, conPolyZono/conZonotope

% Authors:       Niklas Kochdumper
% Written:       06-December-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input arguments are checked in conPolyZono/conZonotope function!

P = polytope(conZonotope(cPZ,varargin{:}));

% set properties
P.bounded.val = true;

% ------------------------------ END OF CODE ------------------------------
