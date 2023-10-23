function han = plot(cPZ,varargin)
% plot - plots an over-approximative projection of a constrained polynomial
%    zonotope
%
% Syntax:
%    han = plot(cPZ)
%    han = plot(cPZ,dims)
%    han = plot(cPZ,dims,type)
%
% Inputs:
%    cPZ - conPolyZono object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'Splits',splits> - number of splits for refinement (default: 8)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    E = [1 0 3;0 1 1;0 0 0];
%    A = [1 -2 1 -2 -0.5];
%    b = 3/2;
%    EC = [2 1 0 0 0; 0 0 2 1 0; 0 0 0 0 1];
%    GI = [0;0.5];
%    
%    cPZ = conPolyZono(c,G,E,A,b,EC,GI);
%     
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Splits',8);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/plot

% Authors:       Niklas Kochdumper
% Written:       19-January-2020
% Last update:   25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[cPZ,dims,NVpairs,splits] = aux_parseInput(cPZ,varargin{:});

% 2. preprocess for plotting
[cPZ,dims] = aux_preprocess(cPZ,dims);

% 3. plot n-dimensional set
han = aux_plotNd(cPZ,dims,NVpairs,splits);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [cPZ,dims,NVpairs,splits] = aux_parseInput(cPZ,varargin)
    % parse input 

    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);

    % check input args
    inputArgsCheck({{cPZ, 'att', 'conPolyZono'},
        {dims,'att','numeric',{'nonempty','vector','integer','positive'}}})
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % read additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end));
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',8);
end

function [cPZ,dims] = aux_preprocess(cPZ,dims)
    % preprocess for plotting 

    % project to desired dimension
    cPZ = project(cPZ,dims);
    dims = 1:length(dims);
end

function han = aux_plotNd(cPZ,dims,NVpairs,splits)
    % plot n-dimensional set

    if length(dims) == 1 % 1d
        % compute enclosing interval
        han = plot(interval(cPZ),dims,NVpairs{:});
    
    elseif length(dims) == 2 % 2d
        % compute enclosing polygon
        pgon = polygon(cPZ,splits);
    
        % plot the polygon
        han = plot(pgon,dims,NVpairs{:});
    
    else % 3d
        han = aux_plot3d(cPZ,dims,NVpairs,splits);
    end
end

function han = aux_plot3d(cPZ,dims,NVpairs,splits)
    % transform to equivalent higher-dimensional polynomial zonotope
    c = [cPZ.c; -cPZ.b];
    G = blkdiag(cPZ.G,cPZ.A);
    E = [cPZ.E,cPZ.EC];
    pZ = polyZonotope(c,G,[],E);
    
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
    % correct dims
    dims = 1:3;

    % plot sets
    han = plotMultipleSetsAsOne(Zs,dims,NVpairs);
end

function res = aux_intersectsNullSpace(obj)
% test if the split set violates the constraints (if it not intersects any
% of the hyperplanes)

    res = true;
    n = length(obj.c);

    % loop over all constraint dimensions
    for i = 4:n
       
        c = zeros(n,1); c(i) = 1;
        hyp = conHyperplane(c,0);
        
        if ~isIntersecting_(hyp,zonotope(obj),'exact')
           res = false;
           return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
