function han = plotZono(cZ,varargin)
% plotZono - Visualizes a 2D-projection of the constraint zonotope and the
%    zonotope without constraints
%
% Syntax:
%    plotZono(cZ)
%    plotZono(cZ,dims,plotOptZ,plotOptCon)
%
% Inputs:
%    cZ - conZonotope object
%    dims - (optional) dimensions of the projection
%    plotOptZ - (optional) cell-array containing the plot settings
%                 for the original zonotope
%    plotOptCon - (optional) cell-array containing the plot settings
%                     for the constrained zonotope
%
% Outputs:
%    han - handle of graphics object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    plotOptZ = {'r','LineWidth',2};
%    plotOptCon = {'b','EdgeColor','none'};
%    plotZono(cZ,[1,2],plotOptZ,plotOptCon);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[cZ,dims,plotOptZ,plotOptCon] = aux_parseInput(cZ,varargin{:});

% 2. preprocess
Z = aux_preprocess(cZ);

% 3. plot n-dimensional set
han = aux_plotNd(cZ,Z,dims,plotOptZ,plotOptCon);

% 4. clean han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [cZ,dims,plotOptZ,plotOptCon] = aux_parseInput(cZ,varargin)
    % parse input arguments
    [dims,plotOptZ,plotOptCon] = setDefaultValues(...
        {[1,2],'b',{'FaceColor','r'}},varargin);
    
    % check input arguments
    inputArgsCheck({{cZ,'att','conZonotope'};
                    {dims,'att','numeric','nonnan'};
                    {plotOptZ,'att','cell'};
                    {plotOptCon,'att','cell'}});
end

function Z = aux_preprocess(cZ)
    % preprocess
    Z = zonotope(cZ.c,cZ.G);
end

function han = aux_plotNd(cZ,Z,dims,plotOptZ,plotOptCon)
    % get hold status 
    holdStatus = ishold;
    
    % plot the original zonotope
    plot(Z,dims,plotOptZ{:});
    
    % set hold to on for plot with constraints
    hold on;
    
    % plot the constrained zonotope and return handle
    han = plot(cZ,dims,plotOptCon{:});
    
    % restore hold status
    if holdStatus
        hold on;
    else
        hold off;
    end
end

% ------------------------------ END OF CODE ------------------------------
