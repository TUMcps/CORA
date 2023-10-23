function cZsplit = split(cZ,varargin)
% split - splits a constrained zonotope into two or more constrained
%    zonotopes that enclose the original conZonotope object
%
% Syntax:
%    cZsplit = split(cZ)
%    cZsplit = split(cZ,splitDim)
%
% Description:
%    If only one input is provided, all possible splits of the
%    interval that over-approximate the constrained zonotope are
%    calculated. If the index of one specific dimension is passed as a
%    second input argument, then the split at this specific dimension is
%    calculated.
%
% Inputs:
%    cZ - conZonotope object
%    splitDim - dimension of the over-approximating interval that is splitted
%
% Outputs:
%    cZsplit - cell array of intervals represented as conZonotope objects
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    cZsplit = split(cZono);
%    
%    figure; xlim([-3,1]); ylim([-3,4]);
%    plot(cZono,[1,2],'FaceColor','r');
%
%    figure; hold on; xlim([-3,1]); ylim([-3,4]);
%    plot(cZsplit{1}{1},[1,2],'FaceColor','b');
%    plot(cZsplit{1}{2},[1,2],'FaceColor','g');
%
%    figure; hold on; xlim([-3,1]); ylim([-3,4]);
%    plot(cZsplit{2}{1},[1,2],'FaceColor','y');
%    plot(cZsplit{2}{2},[1,2],'FaceColor','c');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       27-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% calculate an interval that over-approximates the constrained zonotope
I = interval(cZ);


% Case 1: calculate all possible splits
if nargin == 1
    
    % split all parallelotope generators
    n = length(I);
    cZsplit = cell(n,1);    
    
    for splitDim = 1:n
        cZsplit{splitDim} = aux_splitOneDim(cZ,I,splitDim); 
    end
    
% Case 2: split at the specified generator
elseif nargin == 2
    
    splitDim = varargin{1};
    cZsplit = aux_splitOneDim(cZ,I,splitDim);
    
end

end


% Auxiliary functions -----------------------------------------------------

function cZsplit = aux_splitOneDim(cZ,I,splitDim)

    % interval center and radius
    c = center(I);
    r = rad(I);
    
    % first half of the over-approximating interval
    cTemp = c;
    rTemp = r;
    cTemp(splitDim) = c(splitDim) - 0.5*r(splitDim);
    rTemp(splitDim) = 0.5 * r(splitDim); 
    
    I1 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % second conZonotope object
    cTemp(splitDim) = c(splitDim) + 0.5*r(splitDim);
    I2 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % intersect the intervals with the original conZonotope object
    cZsplit{1} = and_(cZ,I1,'exact');
    cZsplit{2} = and_(cZ,I2,'exact');
    
end
    
% ------------------------------ END OF CODE ------------------------------
