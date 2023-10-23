function useCORAcolors(identifier, varargin)
% useCORAcolors - sets the CORA plotting colormap
%
% Syntax:
%    useCORAcolors(identifier)
%
% Inputs:
%    identifier - name of CORA colormap, one of:
%       - 'CORA:default'             - matlab default color order
%       - 'CORA:contDynamics'        - plot reachSet, initialSet, simRes 
%       - 'CORA:manual'              - CORA manual: default colors
%       - 'CORA:manual-result'       - CORA manual: colors for result plots
%    varargin - depends on identifier, see below
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAcolor, colororder

% Authors:       Tobias Ladner
% Written:       01-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 1
    throw(CORAerror("CORA:notEnoughInputArgs", 1))
end
inputArgsCheck({{identifier, 'str', ...
    {'CORA:default','CORA:contDynamics','CORA:manual','CORA:manual-result'}}})

colors = [];
switch identifier
    case 'CORA:default'
        colors = [ % matlab default colors
            0         0.4470    0.7410 % blue
            0.8500    0.3250    0.0980 % red
            0.9290    0.6940    0.1250 % yellow
            0.4940    0.1840    0.5560 % purple
            0.4660    0.6740    0.1880 % green
            0.3010    0.7450    0.9330 % light blue
            0.6350    0.0780    0.1840 % dark red
        ];
    case 'CORA:contDynamics'
        % varargin - {numColors}
        %    - numColors: number of reachSet colors
        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs', 2))
        end

        numColors = setDefaultValues({1}, varargin);

        reachSetColors = [];
        for cidx = 1:numColors
            reachSetColors = [
                reachSetColors;
                CORAcolor("CORA:reachSet", numColors, cidx)
            ];
        end

        colors = [
            reachSetColors;
            CORAcolor('CORA:initialSet');
            CORAcolor('CORA:simulations');
        ];
    case 'CORA:manual'
        colors = [ % matlab default colors
            0         0.4470    0.7410 % blue
            0.8500    0.3250    0.0980 % red
            0.9290    0.6940    0.1250 % yellow
            0.4940    0.1840    0.5560 % purple
            0.4660    0.6740    0.1880 % green
            0.3010    0.7450    0.9330 % light blue
            0.6350    0.0780    0.1840 % dark red
        ];
    case 'CORA:manual-result'
        colors = [ 
            0.4660    0.6740    0.1880 % green
        ];
end

set(gca(), 'ColorOrder', colors); 
set(gca(), 'ColorOrderIndex', 1);

end

% ------------------------------ END OF CODE ------------------------------
