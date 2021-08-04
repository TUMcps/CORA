function han = plot(hyp,varargin)
% plot - Plots 2-dimensional projection of a constrained hyperplane
%
% Syntax:  
%    plot(hyp)
%    plot(hyp,dims,type)
%
% Inputs:
%    hyp - halfspace object
%    dimensions - (optional) dimensions that should be projected;
%                 assume that other entries of the normal vector are zeros
%    type - (optional) plot settings (LineSpec and name-value pairs)
%    
%
% Outputs:
%    han - handle to the plotted object
%
% Example: 
%    C = [-1 -1; 1 0;-1 0; 0 1; 0 -1];
%    d = [2;3;2;3;2];
%    hyp = conHyperplane(halfspace([1;1],2),C,d);
% 
%    figure; hold on;
%    plot(mptPolytope(C,d),[1,2],'g');
%    xlim([-4,4]);
%    ylim([-4,4]);
%    plot(hyp,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/plot

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    dims = [1,2];
    type{1} = 'b';

    if nargin >= 2 && ~isempty(varargin{1})
        dims = varargin{1};
    end

    if nargin >= 3 && ~isempty(varargin{2})
        type = varargin(2:end); 
    end
    
    % check dimension
    if length(dims) < 2
        error('At least 2 dimensions have to be specified!');
    elseif length(dims) > 3
        error('Only up to 3 dimensions can be plotted!');
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

    if length(dims) == 2
        C = [C;eye(2);-eye(2)];
        d = [d;xLim(2);yLim(2);-xLim(1);-yLim(1)];
    else
        zLim = get(gca,'Zlim');
        C = [C;eye(3);-eye(3)];
        d = [d;xLim(2);yLim(2);zLim(2);-xLim(1);-yLim(1);-zLim(1)];
    end

    vert = lcon2vert(C,d,hyp.h.c(dims)',hyp.h.d);

    % plot constrained hyperplane
    if length(dims) == 2
        han = plot([vert(1,1),vert(2,1)],[vert(1,2),vert(2,2)],type{:});
    else
        % state space transformation
        B = gramSchmidt(hyp.h.c);
        vert_ = B'*vert';
        ind = convhull(vert_(2:end,:)');
        
        han = fill3(vert(ind,1),vert(ind,2),vert(ind,3),type{:}); 
    end

%------------- END OF CODE --------------