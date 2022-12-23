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

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      13-March-2019
% Last update:  14-July-2020 (merge with plotFilled)
%               12-March-2021
%               19-May-2022 (plot marker if ellipsoid is point)
%               25-May-2022 (TL: 1D Plotting)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
dims = setDefaultValues({[1,2]},varargin);

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {dims,'att','numeric',{'nonnan','vector','integer','positive'}}});

% read out plot options
NVpairs = readPlotOptions(varargin(2:end));

% readout 'FaceColor' to decide plot/fill call where necessary
[~,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% only center remaining
if rank(E)==0
    [~,color] = readNameValuePair(NVpairs,'Color',{},'k');
    [~,marker] = readNameValuePair(NVpairs,'Marker',{},'x');
    [~,markerSize] = readNameValuePair(NVpairs,'MarkerSize',{},6);
    if length(dims) == 1
        han = plot(E.q(dims(1)),0,'Color',color,...
                   'Marker',marker,'MarkerSize',markerSize);
    elseif length(dims) == 2
        han = plot(E.q(dims(1)),E.q(dims(2)),'Color',color,...
                   'Marker',marker,'MarkerSize',markerSize);
    elseif length(dims) == 3
        han = plot(E.q(dims(1)),E.q(dims(2),E.q(dims(3))),'Color',color,...
                   'Marker',marker,'MarkerSize',markerSize);
    end
    return;
end

% 1D, 2D vs 3D plot
if length(dims) <= 2
    N = 1000;
    % project ellipsoid
    E_p = project(E,dims);
    % compute points on boundary
    Y = boundary(E_p,N);
    % repeat first entry
    Y = [Y,Y(:,1)];
    %plot and output the handle
    if length(dims) == 1
       han = plot(Y(1,:),zeros(size(Y(1,:))),NVpairs{:});
    else
        if isempty(facecolor) || strcmp(facecolor,'none')
            han = plot(Y(1,:),Y(2,:),NVpairs{:});
        else
            han = fill(Y(1,:),Y(2,:),facecolor,NVpairs{:});
        end
    end

else
    % enclose ellipsoid with zonotope
    E = project(E,dims);
    Z = zonotope(E,100,'outer:norm');
    
    % plot zonotope enclosure
    han = plot(Z,[1,2,3],NVpairs{:});
end

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------