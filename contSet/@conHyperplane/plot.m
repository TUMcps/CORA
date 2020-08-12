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

    C = [C;eye(2);-eye(2)];
    d = [d;xLim(2);yLim(2);-xLim(1);-yLim(1)];

    vert = lcon2vert(C,d,hyp.h.c(dims)',hyp.h.d);

    % plot constrained hyperplane
    han = plot([vert(1,1),vert(2,1)],[vert(1,2),vert(2,2)],type{:});

%------------- END OF CODE --------------