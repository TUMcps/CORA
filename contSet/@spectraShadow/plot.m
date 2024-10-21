function han = plot(SpS,varargin)
% plot - plots a projection of a spectrahedral shadow
%
% Syntax:
%    spec = plot(SpS)
%    spec = plot(SpS,dims)
%    spec = plot(SpS,dims,type)
%
% Inputs:
%    SpS - spectraShadow object 
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs),
%        including added pairs:
%          <'Splits',splits> - number of splits for refinement
%          <'ApproxType',approxType> - 'inner' (default) yields an
%          under-approximative plot, while 'outer' yields an
%          over-approximative plot of the spectrahedral shadow
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    plot(SpS)
%
% Other m-files required: spectraShadow.m
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg, Tobias Ladner
% Written:       02-August-2023
% Last update:   15-October-2024 (TL, extracted to spectraShadow/polytope + simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[SpS,dims,NVpairs,splits,approxType] = aux_parseInput(SpS,varargin{:});

% 2. project
SpS = project(SpS,dims);

% 3. convert to polytope
P = polytope(SpS,approxType,splits);

% 4. plot the polytope
han = plot(P,1:dim(P),NVpairs{:});

% 5. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [SpS,dims,NVpairs,splits,approxType] = aux_parseInput(SpS,varargin)
    % parse input

    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);
    inputArgsCheck({{SpS,'att','spectraShadow'};
        {dims,'att','numeric','vector'}
    });
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % read additional name-value pairs
    NVpairs = varargin(2:end);
    % read out 'Splits', default value given
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',100);
    % read out 'approxType', default value given
    [NVpairs,approxType] = readNameValuePair(NVpairs,'ApproxType','ischar','inner');
    
    % values for name-value pairs are checked in spectraShadow/polytope

end

% ------------------------------ END OF CODE ------------------------------
