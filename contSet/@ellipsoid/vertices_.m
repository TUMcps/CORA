function V = vertices_(E,varargin)
% vertices_ - computes the vertices of a ellipsoid
%
% Syntax:
%    V = vertices_(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    V - vertices
%
% Example: 
%    E = ellipsoid(1,1);
%    V = vertices(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Mark Wetzlinger
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only support 1D ellipsoids
if dim(E) > 1
    throw(CORAerror('CORA:notSupported',...
        'Vertex computation of ellipsoids only supported for 1D cases.'));
end

% compute minimum and maximum using support function
[~,Vmin] = supportFunc_(E,-1,'upper');
[~,Vmax] = supportFunc_(E,1,'upper');

if withinTol(Vmin,Vmax)
    % only pick one vertex
    V = Vmin;
else
    V = [Vmin,Vmax];
end

% ------------------------------ END OF CODE ------------------------------
