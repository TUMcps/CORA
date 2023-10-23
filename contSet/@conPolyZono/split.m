function res = split(cPZ,varargin)
% split - splits a constrained polynomial zonotope along one dimension
%
% Description:
%    If only one input is provided, all possible splits are calculated
%    If the index of one specific dimension is passed as a
%    second input argument, then the split at this specific dimension is
%    calculated.
%
% Syntax:
%    res = split(cPZ)
%    res = split(cPZ, splitDim)
%
% Inputs:
%    cPZ - conPolyZono object
%    splitDim - dimension of the dimension that is splitted
%
% Outputs:
%    res - cell-array containing the splitted objects
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 2 1 2];
%    E = [1 2 1 0; 0 0 1 2; 0 0 0 0];
%    A = [1 1 0.5];
%    b = 0.5;
%    EC = [0 1 0;1 0 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    res = split(cPZ,1);
%
%    figure; hold on;
%    plot(res{1},[1,2],'FaceColor','b','Splits',12);
%    plot(res{2},[1,2],'FaceColor','g','Splits',12);
%    plot(cPZ,[1,2],'r','Splits',12,'LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/split

% Authors:       Niklas Kochdumper
% Written:       08-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% calculate an interval that encloses the constrained polynomial zonotope
int = interval(cPZ);

% Case 1: calculate all possible splits
if nargin == 1
   
    n = length(int);
    res = cell(n,1);    
    
    for splitDim = 1:n
        res{splitDim} = aux_splitOneDim(cPZ,int,splitDim); 
    end
    
% Case 2: split at the specified dimension
elseif nargin == 2
    
    splitDim = varargin{1};
    res = aux_splitOneDim(cPZ,int,splitDim);
    
end

end


% Auxiliary functions -----------------------------------------------------

function cZsplit = aux_splitOneDim(cZ,inter,splitDim)

    % interval center and radius
    c = center(inter);
    r = rad(inter);
    
    % first half of the enclosing interval
    cTemp = c;
    rTemp = r;
    cTemp(splitDim) = c(splitDim) - 0.5*r(splitDim);
    rTemp(splitDim) = 0.5 * r(splitDim); 
    
    I1 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % second half of the enclosing interval 
    cTemp(splitDim) = c(splitDim) + 0.5*r(splitDim);
    I2 = interval(cTemp-rTemp,cTemp+rTemp);
    
    % intersect the intervals with the original conPolyZono object
    cZsplit{1} = and_(cZ,I1,'exact');
    cZsplit{2} = and_(cZ,I2,'exact');
    
end
    
% ------------------------------ END OF CODE ------------------------------
