function [cZsplit] = split(obj,varargin)
% split - splits a constrained zonotope into two or more constrained
%         zonotopes that enclose the original conZonotope object
%
% Syntax:  
%    [cZsplit] = split(obj)
%    [cZsplit] = split(obj, dim)
%
% Description:
%    If only one input is provided, all possible splits of the
%    interval that over-approximates the constrained zonotope are
%    calculated. If the index of one specific dimension is passed as a
%    second input argument, then the split at this specific dimension is
%    calculated.
%
% Inputs:
%    obj - conZonotope object
%    dim - dimension of the over-approximating interval that is splitted
%
% Outputs:
%    pZsplit - cell array of intervals represented as conZonotope objects
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    cZsplit = split(cZono);
%    
%    figure
%    plotFilled(cZono,[1,2],'r');
%    xlim([-3,1]);
%    ylim([-3,4]);
%
%    figure
%    hold on
%    plotFilled(cZsplit{1}{1},[1,2],'b');
%    plotFilled(cZsplit{1}{2},[1,2],'g');
%    xlim([-3,1]);
%    ylim([-3,4]);
%
%    figure
%    hold on
%    plotFilled(cZsplit{2}{1},[1,2],'y');
%    plotFilled(cZsplit{2}{2},[1,2],'c');
%    xlim([-3,1]);
%    ylim([-3,4]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      27-July-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

% calculate an interval that over-approximates the constrained zonotope
inter = interval(obj);


% Case 1: calculate all possible splits
if nargin == 1
    
    % split all parallelotope generators
    n = length(inter);
    cZsplit = cell(n,1);    
    
    for dim = 1:n
        cZsplit{dim} = splitOneDim(obj,inter,dim); 
    end
    
% Case 2: split at the specified generator
elseif nargin == 2
    
    dim = varargin{1};
    cZsplit = splitOneDim(obj,inter,dim);
    
end



% Auxiliary functions -----------------------------------------------------

function cZsplit = splitOneDim(cZ,inter,dim)

    % interval center and radius
    r = rad(inter);
    c = mid(inter);
    
    % first halp of the over-approximating interval
    cTemp = c;
    rTemp = r;
    cTemp(dim) = c(dim) - 0.5*r(dim);
    rTemp(dim) = 0.5 * r(dim); 
    
    int1 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % second conZonotope object
    cTemp(dim) = c(dim) + 0.5*r(dim);
    int2 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % intersect the intervals with the original conZonotope object
    cZsplit{1} = cZ & int1;
    cZsplit{2} = cZ & int2;
    
    
    
%------------- END OF CODE --------------