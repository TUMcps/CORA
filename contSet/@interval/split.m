function res = split(I,n)
% split - split an interval in one dimension
%
% Syntax:
%    res = split(I,n)
%
% Inputs:
%    I - interval object
%    n - index of the dimension that is splitted
%
% Outputs:
%    res - cell-array containing the splitted intervals
%
% Example: 
%    I = interval([-1;-1],[1;1]);
%    I_split = split(I,1);
% 
%    figure; hold on; xlim([-2,2]); ylim([-2,2]);
%    plot(I_split{1},[1,2],'FaceColor','g');
%    plot(I_split{2},[1,2],'FaceColor','b');
%    plot(I,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/split

% Authors:       Niklas Kochdumper
% Written:       25-July-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

m = center(I);

sup = I.sup;
infi = I.inf;

sup(n) = m(n);
infi(n) = m(n);

res = {interval(I.inf,sup),interval(infi,I.sup)};

% ------------------------------ END OF CODE ------------------------------
