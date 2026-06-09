function han = plotHalfspaces(P,varargin)
% plotHalfspaces - plots all halfspaces of a polytope
%
% Syntax:
%    han = plotHalfspaces(P)
%    han = plotHalfspaces(P,dims)
%    han = plotHalfspaces(P,dims,varargin)
%
% Inputs:
%    P - polytope
%    dims - (optional) dimensions for projection
%    varargin - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    P = polytope.generateRandom('Dimension',2,'IsBounded',true);
%    
%    figure; hold on;
%    plot(P); enlargeAxis;
%    plotHalfspaces(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Tobias Ladner
% Written:       22-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
dims = setDefaultValues({1:2},varargin);
numHalfspaces = numel(P.b);
if numHalfspaces == 0
    CORAwarning('CORA:plot','Polytope does not contain any halfspaces.')
    return
end

% extract all halfspaces
HSs = arrayfun(@(i) polytope(P.A(i,:),P.b(i)), 1:numHalfspaces, 'UniformOutput',false);

% plot halfspaces 
% (face alpha is set s.t P shows up as if plotted with 'Filled', true
han = plotMultipleSetsAsOne(HSs,dims,[{'FaceAlpha',0.2/numHalfspaces},varargin(2:end)]);

% clear output
if nargout == 0
    clear han
end

% ------------------------------ END OF CODE ------------------------------
