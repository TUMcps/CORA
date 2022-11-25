function res = split(obj,dim)
% split - split an interval in one dimension
%
% Syntax:  
%    res = split(obj,dim)
%
% Inputs:
%    obj - interval object
%    dim - index of the dimension that is splitted
%
% Outputs:
%    res - cell-array containing the splitted intervals
%
% Example: 
%    int = interval([-1;-1],[1;1]);
%    intSplit = split(int,1);
% 
%    figure
%    hold on
%    plot(intSplit{1},[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(intSplit{2},[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(int,[1,2],'r');
%    xlim([-2,2]);
%    ylim([-2,2]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/split

% Author:       Niklas Kochdumper
% Written:      25-July-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

m = center(obj);

sup = obj.sup;
infi = obj.inf;

sup(dim) = m(dim);
infi(dim) = m(dim);

res = {interval(obj.inf,sup),interval(infi,obj.sup)};

%------------- END OF CODE --------------