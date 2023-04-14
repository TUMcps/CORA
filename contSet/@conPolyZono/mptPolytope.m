function P = mptPolytope(cPZ,varargin)
% mptPolytope - Encloses a constrained polynomial zonotope by a polytope
%
% Syntax:  
%    P = mptPolytope(cPZ)
%    P = mptPolytope(cPZ,type)
%    P = mptPolytope(cPZ,type,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    type - algorithm used to compute polytope enclosure ('linearize',
%           'extend', or 'all')
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    expMat_ = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    expMat = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    P = mptPolytope(cPZ);
%
%    figure; hold on
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(P,[1,2],'b');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope, conPolyZono/conZonotope

% Author:       Niklas Kochdumper
% Written:      06-December-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% input arguments are checked in conPolyZono/conZonotope function!

P = mptPolytope(conZonotope(cPZ,varargin{:}));

%------------- END OF CODE --------------