function han = plot3D(pZ,varargin)
% plot3D - plots a 3D projection of a polynomial zonotope
%
% Syntax:
%    han = plot3D(pZ)
%    han = plot3D(pZ,NVpairsPlot)
%
% Inputs:
%    pZ - projected polyZonotope object
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
%    NVpairsVertices - (optional) vertex computation settings
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Tobias Ladner
% Written:       29-March-2018
% Last update:   23-June-2020 (MW, harmonize with other plot functions)
%                14-July-2020 (MW, merge with plotFilled)
%                25-May-2022 (TL, 1D Plotting)
%                23-February-2023 (TL, enlarge polygon if 2 regions)
%                23-March-2023 (TL, bugfix Splits=0, clean up)
%                05-April-2023 (TL, clean up using plotPolygon)
%                11-July-2023 (TL, bug fix holes)
%                13-April-2024 (TL, bug fix 3d)
%                17-July-2024 (TL, minor speed up)
% Last revision: 12-July-2023 (TL, restructure)
%                15-October-2024 (TL, split into plot1D/plot2D/plot3D)

% ------------------------------ BEGIN CODE -------------------------------

% parse additional parameters
[NVpairsPlot,NVpairsVertices] = setDefaultValues({{},{}},varargin);

% read splits
splits = NVpairsVertices{1};

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape
pZsplit{1} = pZ;

for i=1:splits
    pZnew = cell(2*length(pZsplit),1);
    for j=1:length(pZsplit)
        res = splitLongestGen(pZsplit{j});
        pZnew{2*j-1} = res{1};
        pZnew{2*j} = res{2};
    end
    pZsplit = pZnew;
end

% convert all sets to a zonotope
Zs = cell(1,numel(pZsplit));
for i = 1:numel(pZsplit)
    Zs{i} = zonotope(pZsplit{i});
end

% plot all sets as one
han = plotMultipleSetsAsOne(Zs,1:3,NVpairsPlot);

end

% ------------------------------ END OF CODE ------------------------------
