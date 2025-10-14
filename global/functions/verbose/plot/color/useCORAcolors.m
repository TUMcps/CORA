function useCORAcolors(identifier, varargin)
% useCORAcolors - sets the CORA plotting colormap
%
% Syntax:
%    useCORAcolors(identifier)
%
% Inputs:
%    identifier - name of CORA colormap, one of:
%       - 'CORA:default'             - matlab default color order
%       - 'CORA:contDynamics'        - plot reachSet, initialSet, traj 
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
% Last update:   27-August-2025 (TL, added matlab default colors)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,Inf);
inputArgsCheck({{identifier, 'str', ...
    {'CORA:default','MATLAB:default','CORA:contDynamics','CORA:manual','CORA:manual-result'}}})

colors = [];
switch identifier
    case 'CORA:default'
        colors = [ % cora default colors (differ from matlab since R2025a)
            CORAcolor('CORA:color1')
            CORAcolor('CORA:color2')
            CORAcolor('CORA:color3')
            CORAcolor('CORA:color4')
            CORAcolor('CORA:color5')
            CORAcolor('CORA:color6')
            CORAcolor('CORA:color7')
        ];
    case 'MATLAB:default'
        colors = [ % matlab default colors
            CORAcolor('MATLAB:color1')
            CORAcolor('MATLAB:color2')
            CORAcolor('MATLAB:color3')
            CORAcolor('MATLAB:color4')
            CORAcolor('MATLAB:color5')
            CORAcolor('MATLAB:color6')
            CORAcolor('MATLAB:color7')
        ];
    case 'CORA:contDynamics'
        % varargin - {numColors}
        %    - numColors: number of reachSet colors
        narginchk(1,2);

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
            CORAcolor('CORA:color1')
            CORAcolor('CORA:color2')
            CORAcolor('CORA:color3')
            CORAcolor('CORA:color4')
            CORAcolor('CORA:color5')
            CORAcolor('CORA:color6')
            CORAcolor('CORA:color7')
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
