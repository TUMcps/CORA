function han = plot(hs,varargin)
% plot - plots a projection of a halfspace
%
% Syntax:
%    han = plot(hs)
%    han = plot(hs,dims)
%    han = plot(hs,dims,type)
%
% Inputs:
%    hs - halfspace object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    hs = halfspace([1;1],0);
% 
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(hs,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/plot

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       23-August-2013
% Last update:   19-November-2019 (NK, plot area instead of line)
%                25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
%                26-July-2023 (TL, getUnboundedAxisLimits)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[hs,dims,NVpairs] = aux_parseInput(hs,varargin{:});

% 2. preprocess (and convert to polytope)
[P,dims] = aux_preprocess(hs,dims);

% 3. plot polytope
han = plot(P,dims,NVpairs{:});

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [hs,dims,NVpairs] = aux_parseInput(hs,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{hs,'att','halfspace'};
                    {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % read additional name-value pairs
    if size(dims) == 1
        NVpairs = readPlotOptions(varargin(2:end),'contour');
    else
        NVpairs = readPlotOptions(varargin(2:end),'fill');
    end
end

function [P,dims] = aux_preprocess(hs,dims)

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    
    % convert to polytope
    % projection with other dimensions = 0
    % and intersection with currently visible axis
    if length(dims) == 1
        C = [hs.c(dims)';1;-1];
        d = [hs.d;xLim(2);-xLim(1)];
    elseif length(dims) == 2
        C = [hs.c(dims)';eye(2);-eye(2)];
        d = [hs.d;xLim(2);yLim(2);-xLim(1);-yLim(1)];
    else
        C = [hs.c(dims)';eye(3);-eye(3)];
        d = [hs.d;xLim(2);yLim(2);zLim(2);-xLim(1);-yLim(1);-zLim(1)];
    end
    P = polytope(C,d);
end

% ------------------------------ END OF CODE ------------------------------
