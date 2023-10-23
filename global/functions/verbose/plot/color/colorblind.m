function rgb = colorblind(color)
% colorblind - colors suitable for colorblind people
%
% Syntax:
%    rgb = colorblind('b')
%
% Inputs:
%    color - desired color, currently only
%               'blue', 'b' for blue,
%               'red', 'r' for red,
%               'yellow', 'y' for yellow, and
%               'gray' for gray supported         
%
% Outputs:
%    rgb - red-green-blue values of the desired color
% 
% References:
%   [1] T.B. Plante, M. Cushman. Choosing color palettes for scientific
%       figures. Research and practice in thrombosis and haemostasis 4(2),
%       pp. 176-180, 2020.

% Authors:       Mark Wetzlinger
% Written:       24-September-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% provided colors
colors = {'blue','b','red','r','yellow','y','green','g','gray'};
if ~any(ismember(colors,color))
    throw(CORAerror('CORA:wrongValue','second','Desired color not provided!'));
end

rgb = [0,0,0];
if any(ismember({'blue','b'},color))
    rgb = [0, 92, 171] ./ 255;
elseif any(ismember({'red','r'},color))
    rgb = [227, 27, 35] ./ 255;
elseif any(ismember({'yellow','y'},color))
    rgb = [255, 195, 37] ./ 255;
elseif any(ismember({'green','g'},color))
    rgb = [0.4660 0.6740 0.1880];
elseif any(ismember({'gray'},color))
    rgb = [230, 241, 238] ./ 255;
end

% ------------------------------ END OF CODE ------------------------------
