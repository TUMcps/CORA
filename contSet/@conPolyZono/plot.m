function han = plot(cPZ,varargin)
% plot - plots enclosure of the 2D or 3D projection of a constrained 
%        polynomial zonotope
%
% Syntax:  
%    han = plot(cPZ)
%    han = plot(cPZ,dims,linespec)
%    han = plot(cPZ,dims,linespec,'Splits',splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    dims - dimensions that should be projected
%    linespec - (optional) LineSpec properties
%    splits - (optional) number of splits for refinement
%    type - (optional) name-value pairs
%
% Outputs:
%    han - handle for the resulting graphics object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    expMat = [1 0 3;0 1 1;0 0 0];
%    A = [1 -2 1 -2 -0.5];
%    b = 3/2;
%    expMat_ = [2 1 0 0 0; 0 0 2 1 0; 0 0 0 0 1];
%    Grest = [0;0.5];
%    
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
%     
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Splits',8);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/plot

% Author:       Niklas Kochdumper
% Written:      19-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % default values for the optional input arguments
    dims = [1,2];
    linespec = 'b';
    splits = 8;
    NVpairs = {};
    filled = 0;

    % parse input arguments
    if nargin > 1 && ~isempty(varargin{1})
        dims = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
        % read additional name-value pairs
        [linespec,NVpairs] = readPlotOptions(varargin(2:end));
        [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical',filled);
        [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',splits);
    end

    % check dimension
    if length(dims) < 2
        error('At least 2 dimensions have to be specified!');
    elseif length(dims) > 3
        error('Only up to 3 dimensions can be plotted!');
    end

    % project to desired dimensions
    cPZ = project(cPZ,dims);

    % 2D vs 3D plot
    if length(dims) == 2

        % compute enclosing polygon
        pgon = polygon(cPZ,splits);

        % plot the polygon
        if filled
            han = plot(pgon,[1,2],linespec,NVpairs{:},'Filled',true);
        else
            han = plot(pgon,[1,2],linespec,NVpairs{:});
        end

    else

        % transform to equivalent higher-dimensional polynomial zonotope
        c = [cPZ.c; -cPZ.b];
        G = blkdiag(cPZ.G,cPZ.A);
        expMat = [cPZ.expMat,cPZ.expMat_];

        pZ = polyZonotope(c,G,[],expMat);
        
        % split the polynomial zonotope multiple times to obtain a better 
        % over-approximation of the real shape
        pZsplit{1} = pZ;

        for i = 1:splits
            pZnew = [];
            for j = 1:length(pZsplit)
                res = splitLongestGen(pZsplit{j});
                if intersectsNullSpace(res{1})
                    pZnew{end+1} = res{1};
                end
                if intersectsNullSpace(res{2})
                    pZnew{end+1} = res{2};
                end
            end
            pZsplit = pZnew;
        end

        % check if set is empty
        if isempty(pZsplit)
           error('Set is empty!'); 
        end
        
        % loop over all parallel sets
        hold on;
        for i = 1:length(pZsplit)
            
            zono = zonotope(project(pZsplit{i},[1,2,3]));
            zono = zonotope([zono.Z,cPZ.Grest]);
            
            if filled
                han = plot(zono,[1,2,3],linespec,NVpairs{:},'Filled',true); 
            else
                han = plot(zono,[1,2,3],linespec,NVpairs{:});  
            end
        end
    end
end

% Auxiliary Functions -----------------------------------------------------

function res = intersectsNullSpace(obj)
% test if the split set violates the constraints (if it not intersects any
% of the hyperplanes)

    res = 1;
    n = length(obj.c);

    % loop over all constraint dimensions
    for i = 4:n
       
        c = zeros(n,1); c(i) = 1;
        hs = conHyperplane(c,0);
        
        if ~isIntersecting(hs,zonotope(obj))
           res = 0;
           return;
        end
    end
end

%------------- END OF CODE --------------