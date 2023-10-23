function han = plot(I,varargin)
% plot - plots a projection of an interval 
%
% Syntax:
%    han = plot(I)
%    han = plot(I,dims)
%    han = plot(I,dims,type)
%
% Inputs:
%    I - interval object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    plot(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       31-July-2016
% Last update:   12-July-2023 (TL, use vertices)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[I,dims,NVpairs] = aux_parseInput(I,varargin{:});

% 2. preprocess
V = aux_preprocess(I,dims);

% 3. plot
han = plotPolygon(V,'ConvHull',true,NVpairs{:});

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [I,dims,NVpairs] = aux_parseInput(I,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{I,'att','interval'};
                    {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % read additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end));
end

function V = aux_preprocess(I,dims)
    % preprocess

    % project
    I = project(I,dims);

    % compute vertices
    V = vertices(I);

    % infinity values in V will be considered in plotPolygon    
end

% ------------------------------ END OF CODE ------------------------------
