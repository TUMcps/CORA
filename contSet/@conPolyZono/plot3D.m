function han = plot3D(cPZ,varargin)
% plot3D - plots a 3D projection of a polynomial zonotope
%
% Syntax:
%    han = plot3D(cPZ)
%    han = plot3D(cPZ,NVpairsPlot)
%
% Inputs:
%    cPZ - projected conPolyZono object
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

% transform to equivalent higher-dimensional polynomial zonotope
c = [cPZ.c; -cPZ.b];
G = blkdiag(cPZ.G,cPZ.A);
E = [cPZ.E,cPZ.EC];
pZ = polyZonotope(c,G,[],E);

% correct dims
dims = 1:3;

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape
pZsplit{1} = pZ;
for i = 1:splits
    pZnew = [];
    for j = 1:length(pZsplit)
        % split pZ
        res = splitLongestGen(pZsplit{j});

        % check contrains of splitted sets
        if aux_intersectsNullSpace(res{1})
            pZnew{end+1} = res{1};
        end
        if aux_intersectsNullSpace(res{2})
            pZnew{end+1} = res{2};
        end
    end
    pZsplit = pZnew;
end

% check if set is empty
if isempty(pZsplit)
     throw(CORAerror('CORA:emptySet'));
end

% convert all sets into a zonotope
Zs = cell(1,length(pZsplit));
for i=1:length(pZsplit)
    % read pZ i
    pZi = pZsplit{i};
    
    % project to correct dims
    pZi = project(pZi,dims);

    % convert to zonotope
    Zi = zonotope(pZi);

    % add independent generators
    Zi = zonotope(Zi.c, [Zi.G,cPZ.GI]);

    % add to list
    Zs{i} = Zi;
end

% plot sets
han = plotMultipleSetsAsOne(Zs,dims,NVpairsPlot);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_intersectsNullSpace(pZ)
% test if the split set violates the constraints (if it not intersects any
% of the hyperplanes)

    res = true;
    n = dim(pZ);

    % loop over all constraint dimensions
    for i = 4:n
        P = polytope([],[],unitvector(i,n)',0);
        if ~isIntersecting_(P,zonotope(pZ),'exact',1e-8)
            res = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
