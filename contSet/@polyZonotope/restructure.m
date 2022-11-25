function res = restructure(pZ,method,order,varargin)
% restructure - Calculate a new over-approxmiating representation of a
%               polynomial zonotope in such a way that there remain no
%               independent generators
%
% Syntax:  
%    res = restructure(pZ, method, order)
%    res = restructure(pZ, method, order, genOrder)
%
% Inputs:
%    pZ - polyZonotope object
%    method - method used to calculate the new representation ('zonotope'
%             or 'reduce', 'reduceFull', or 'reducePart')
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    genOrder - desired zonotope order of the resulting polynomial zonotope
%               (only for method = 'reduce...')
%
% Outputs:
%    res - resuling polyZonotope object which over-approximates pZ
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'zonotopeGirard',2);
%    pZnew2 = restructure(pZ,'reduceGirard',2);
%
%    hold on
%    plot(pZnew2,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    ylim([-10,10]);
%
%    figure
%    hold on
%    plot(pZnew1,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    ylim([-10,10]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reduce

% Author:       Niklas Kochdumper
% Written:      25-July-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    genOrder = inf;
    if nargin >= 4
       genOrder = varargin{1}; 
    end

    % parse string for the method
    if startsWith(method,'zonotope')
        spec = 1;
        redMeth = method(9:end);
    elseif startsWith(method,'reduceFull')
        spec = 2;
        redMeth = method(11:end);
    elseif startsWith(method,'reducePart')
        spec = 3;
        redMeth = method(11:end);
    elseif startsWith(method,'reduce')
        spec = 0;
        redMeth = method(7:end);
    else
        error('Wrong value for input argument method!');
    end
    
    redMeth(1) = lower(redMeth(1));
    
    % restructure with the selected method
    if spec == 1
        res = restructureZono(pZ,order,redMeth);
    elseif spec == 2
        res = restructureReduceFull(pZ,order,redMeth);
    elseif spec == 3
        res = restructureReducePart(pZ,order,redMeth);
    else
        res = restructureReduce(pZ,order,redMeth,genOrder);
    end

%------------- END OF CODE --------------