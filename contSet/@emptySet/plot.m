function han = plot(O,varargin)
% plot - plots a projection of an emptySet
%
% Syntax:
%    han = plot(O)
%    han = plot(O,dims)
%    han = plot(O,dims,type)
%
% Inputs:
%    O - emptySet object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    O = emptySet(2);
% 
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(O,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       03-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[O,dims,NVpairs] = aux_parseInput(O,varargin{:});

% 2. plot
han = aux_plot(O,dims,NVpairs);

% 3. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [O,dims,NVpairs] = aux_parseInput(O,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{O,'att','emptySet'};
                    {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    if max(dims) > dim(O)
        throw(CORAerror("CORA:wrongValue",'second','Specified dimensions must not exceed the dimension of the set.'))
    end
    
    % read additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end),'contour');
end

function han = aux_plot(O,dims,NVpairs)
% plot N-dimensional set

    if length(dims) <= 2
        % plot 1d/2d set
        han = plotPolygon(zeros(length(dims),0),NVpairs{:});
    else
        han = plotPolytope3D(zeros(length(dims),0),NVpairs{:});
    end
    
end

% ------------------------------ END OF CODE ------------------------------
