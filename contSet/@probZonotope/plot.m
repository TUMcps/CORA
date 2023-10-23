function han = plot(probZ,varargin)
% plot - plots a projection of a probabilistic zonotope
%
% Syntax:
%    han = plot(probZ)
%    han = plot(probZ,dims)
%    han = plot(probZ,dims,type)
%
% Inputs:
%    probZ - probZonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'m',m> - m-sigma value (default: probZ.gamma)
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    plot(probZ,[1,2],'FaceColor','red');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       03-August-2007
% Last update:   17-July-2020
%                25-May-2022 (TL, 1D Plotting)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[probZ,dims,NVpairs,m,facecolor] = aux_parseInput(probZ,varargin{:});

% 2. preprocess
eP = aux_preprocess(probZ,m,dims);

% 3. plot n-dimensional set
han = aux_plotNd(eP,NVpairs,facecolor);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [probZ,dims,NVpairs,m,facecolor] = aux_parseInput(probZ,varargin)
    % default values
    dims = setDefaultValues({[1,2]},varargin);
    m = probZ.gamma;
    
    % check input arguments
    inputArgsCheck({{probZ,'att','probZonotope'};
                    {dims,'att','numeric',{'nonempty','integer','positive','vector'}}});
    
    % parse plot options
    NVpairs = readPlotOptions(varargin(2:end),'surf');
    [NVpairs,m] = readNameValuePair(NVpairs,'m','isscalar',m);
    % readout 'FaceColor' to decide plot/fill call where necessary
    [~,facecolor] = readNameValuePair(NVpairs,'FaceColor');
end

function eP = aux_preprocess(probZ,m,dims)
    % preprocess
    
    % one-dimensional case
    if length(dims) == 1
        probZ = project(probZ, dims);
        probZ = [1;0] * probZ;
        dims = [1,2];
    end

    % compute enclosing probability
    eP = enclosingProbability(probZ,m,dims);
end

function han = aux_plotNd(eP,NVpairs,facecolor)
    % plot and output the handle
    if isempty(facecolor) || strcmp(facecolor,'none')
        han = mesh(eP.X,eP.Y,eP.P,NVpairs{:});
    else
        han = surf(eP.X,eP.Y,eP.P,NVpairs{:});
    end
end

% ------------------------------ END OF CODE ------------------------------
