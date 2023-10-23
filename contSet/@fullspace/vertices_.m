function V = vertices_(fs,varargin)
% vertices_ - returns the vertices of a full-dimensional space
%
% Syntax:
%    V = vertices_(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    V - vertices
%
% Example: 
%    fs = fullspace(2);
%    V = vertices(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported',...
        'Vertices computation of R^0 not supported.'));
end

% convert to interval and compute vertices
V = vertices_(interval(fs));

% ------------------------------ END OF CODE ------------------------------
