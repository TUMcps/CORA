function [val,x] = supportFunc_(E,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of an ellipsoid along a
%    certain direction (see Def. 2.1.2 in [1]) 
%
% Syntax:
%    val = supportFunc_(E,dir)
%    [val,x] = supportFunc_(E,dir,type)
%
% Inputs:
%    E - ellipsoid object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower',upper','range')
%
% Outputs:
%    val - bound of the ellipsoid in the specified direction
%    x - point for which holds: dir'*x=val
%
% Example: 
%    E = ellipsoid([5 7;7 13],[1;2]);
%    dir = [1;1];
%
%    [val,x] = supportFunc(E,dir);
%   
%    figure; hold on; box on;
%    plot(E,[1,2],'b');
%    plot(conHyperplane(dir,val),[1,2],'g');
%    plot(x(1),x(2),'.r','MarkerSize',20);
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, zonotope/supportFunc_

% Authors:       Victor Gassmann
% Written:       20-November-2019
% Last update:   12-March-2021
%                27-July-2021 (fixed degenerate case)
%                04-July-2022 (VG, class array case)
%                29-March-2023 (VG, changed to explicit comp of x vector)
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

% set is just a point
if representsa_(E+(-E.q),'origin',eps)
    val = dir'*E.q; x = E.q; return
end

if strcmp(type,'upper')
    val = dir'*E.q + 1*sqrt(dir'*E.Q*dir);
    x = E.q +E.Q*dir/sqrt(dir'*E.Q*dir);
elseif strcmp(type,'lower')
    val = dir'*E.q + (-1)*sqrt(dir'*E.Q*dir);
    x = E.q -E.Q*dir/sqrt(dir'*E.Q*dir);
elseif strcmp(type,'range')
    val = interval(dir'*E.q + (-1)*sqrt(dir'*E.Q*dir),...
        dir'*E.q + 1*sqrt(dir'*E.Q*dir));
    x = [E.q-E.Q*dir/sqrt(dir'*E.Q*dir),...
         E.q+E.Q*dir/sqrt(dir'*E.Q*dir)];
end

% ------------------------------ END OF CODE ------------------------------
