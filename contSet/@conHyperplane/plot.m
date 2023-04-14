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
%    plot(mptPolytope(C,d),[1,2],'g');
%    plot(hyp,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/plot

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  25-May-2022 (TL: 1D Plotting)
%               05-April-2023 (TL: clean up using plotPolygon)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
dims = setDefaultValues({[1,2]},varargin);

% check input arguments
inputArgsCheck({{hyp,'att','conHyperplane'};
                {dims,'att','numeric',{'vector','positive','integer'}}});

% parse plot options
NVpairs = readPlotOptions(varargin(2:end));

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% get size of current plot
xLim = get(gca,'Xlim');
yLim = get(gca,'Ylim');

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
    zLim = get(gca,'Zlim');
    C = [C;eye(3);-eye(3)];
    d = [d;xLim(2);yLim(2);zLim(2);-xLim(1);-yLim(1);-zLim(1)];
end

vert = lcon2vert(C,d,hyp.h.c(dims)',hyp.h.d)';

% plot constrained hyperplane
if length(dims) == 1
    han = plotPolygon(vert,NVpairs{:});
    
elseif length(dims) == 2 
    han = plotPolygon(vert,NVpairs{:});
else
    % state space transformation
    B = gramSchmidt(hyp.h.c);
    vert_ = B'*vert;
    ind = convhull(vert_(2:end,:)');
    
    han = plotPolygon([vert(1,ind);vert(2,ind);vert(3,ind)],NVpairs{:}); 
end

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------