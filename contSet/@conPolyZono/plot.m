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
% Last update:  25-May-2022 (TL: 1D Plotting)
% Last revision:---

%------------- BEGIN CODE --------------

    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);

    % read additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end));
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',8);

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end

    % project to desired dimensions
    cPZ = project(cPZ,dims);

    if length(dims) == 1
        % add zeros to 2nd dimension
        cPZ = cartProd_(cPZ,0,'exact');
        dims = [1;2];
    end

    if length(dims) == 2

        % compute enclosing polygon
        pgon = polygon(cPZ,splits);

        % plot the polygon
        han = plot(pgon,[1,2],NVpairs{:});

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
        
        % loop over all parallel sets
        hold on;
        for i = 1:length(pZsplit)
            
            zono = zonotope(project(pZsplit{i},[1,2,3]));
            zono = zonotope([zono.Z,cPZ.Grest]);
            
            han = plot(zono,[1,2,3],NVpairs{:});  
        end
    end
    
    if nargout == 0
        clear han;
    end
end

% Auxiliary Functions -----------------------------------------------------

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

%------------- END OF CODE --------------