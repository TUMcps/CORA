function han = plot(hyp,varargin)
% plot - plots a projection of a constrained hyperplane
%
% Syntax:
%    han = plot(hyp)
%    han = plot(hyp,dims)
%    han = plot(hyp,dims,type)
%
% Inputs:
%    hyp - conHyperplane object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%    
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    C = [-1 -1; 1 0;-1 0; 0 1; 0 -1];
%    d = [2;3;2;3;2];
%    hyp = conHyperplane(halfspace([1;1],2),C,d);
% 
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(polytope(C,d),[1,2],'g');
%    plot(hyp,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/plot

% Authors:       Niklas Kochdumper
% Written:       19-November-2019
% Last update:   25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
%                26-July-2023 (TL, getUnboundedAxisLimits)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[hyp, dims, NVpairs] = aux_parseInput(hyp,varargin{:});

% 2. preprocess for plotting
[V,dims] = aux_preprocess(hyp,dims);

% 3. plot n-dimensional set
han = aux_plotNd(hyp,V,dims,NVpairs);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [hyp, dims, NVpairs] = aux_parseInput(hyp,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{hyp,'att','conHyperplane'};
                    {dims,'att','numeric',{'vector','positive','integer'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % parse plot options
    NVpairs = readPlotOptions(varargin(2:end));
end

function [V,dims] = aux_preprocess(hyp,dims)
    % preprocess for plotting

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    
    % compute vertices
    if ~isempty(hyp.C)
        C = hyp.C(:,dims);
        d = hyp.d;
    else
        C = []; d = [];
    end
    
    if length(dims) == 1
        C = [C;1;-1];
        d = [d;xLim(2);-xLim(1)];
    elseif length(dims) == 2
        C = [C;eye(2);-eye(2)];
        d = [d;xLim(2);yLim(2);-xLim(1);-yLim(1)];
    else
        C = [C;eye(3);-eye(3)];
        d = [d;xLim(2);yLim(2);zLim(2);-xLim(1);-yLim(1);-zLim(1)];
    end
    
    V = lcon2vert(C,d,hyp.a(dims),hyp.b)';
    dims = 1:length(dims);
end

function han = aux_plotNd(hyp,V,dims,NVpairs)
    % plot n-dimensional set

    if length(dims) <= 2 % 1d, 2d
        han = plotPolygon(V,NVpairs{:});
        
    else % 3d
        % state space transformation
        B = gramSchmidt(hyp.a');
        vert_ = B'*V;
        ind = convhull(vert_(2:end,:)');
        
        han = plotPolygon([V(1,ind);V(2,ind);V(3,ind)],NVpairs{:}); 
    end
end

% ------------------------------ END OF CODE ------------------------------
