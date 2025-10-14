function V = readVerticesFromFigure(obj)
% readVerticesFromFigure - reads out the vertices/points of a given
%    graphics object handle
%
% Syntax:
%    V = readVerticesFromFigure(obj)
%
% Inputs:
%    obj - graphics object handle
%
% Outputs:
%    V - numeric, vertices/points (n x p)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       15-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% as Matlab is a great language, it stores it differently (transposed)
% for different graphics objects 
% -> returning it as all vertices/points are specified in CORA (n x p)
switch class(obj)
    case 'matlab.graphics.primitive.Patch'
        % requires transpose
        V = [obj.XData';obj.YData'];
        if ~isempty(obj.ZData)
            V = [V;obj.ZData'];
        end

    case 'matlab.graphics.chart.primitive.Line'
        % already as in CORA
        V = [obj.XData;obj.YData];
        if ~isempty(obj.ZData)
            V = [V;obj.ZData];
        end

    otherwise
        throw(CORAerror('CORA:specialError',sprintf('Unknown graphics object handle ''%s''', class(obj))))

end

% ------------------------------ END OF CODE ------------------------------
