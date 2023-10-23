function han = plot(E,varargin)
% plot - plots a projection of an ellipsoid
%
% Syntax:
%    han = plot(E)
%    han = plot(E,dims)
%    han = plot(E,dims,type)
%
% Inputs:
%    E - ellipsoid object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    E = ellipsoid([1 0 0; 0 1 0;0 0 3]);
%    plot(E)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       13-March-2019
% Last update:   14-July-2020 (merge with plotFilled)
%                12-March-2021
%                19-May-2022 (plot marker if ellipsoid is point)
%                25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[E,dims,NVpairs] = aux_parseInput(E,varargin{:});

% 2. plot n-dimensional set
han = aux_plotNd(E,dims,NVpairs);

% 3. clearn han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [E,dims,NVpairs] = aux_parseInput(E,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{E,'att','ellipsoid','scalar'};
                    {dims,'att','numeric',{'nonnan','vector','integer','positive'}}});
    
    % read additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end));
end

function han = aux_plotNd(E,dims,NVpairs)
    if rank(E)==0  % only center remaining
        % read center point in desired dimensions
        V = E.q(dims);

        % plot single point
        han = plotPolygon(V,NVpairs{:});
    
    elseif length(dims) <= 2 % 1d, 2d
        N = 1000;
        % project ellipsoid
        E_p = project(E,dims);
        % compute points on boundary
        V = boundary(E_p,N);
        % repeat first entry
        V = [V,V(:,1)];
        % plot and output the handle
        han = plotPolygon(V, NVpairs{:});
    
    else % 3d
        % enclose ellipsoid with zonotope
        E = project(E,dims);
        Z = zonotope(E,100,'outer:norm');
        
        % plot zonotope enclosure
        han = plot(Z,[1,2,3],NVpairs{:});
    end
end

% ------------------------------ END OF CODE ------------------------------
