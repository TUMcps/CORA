function val = supportFunc(E,dir,varargin)
% supportFunc - Calculate the upper or lower bound of an ellipsoid object
%               along a certain direction (see Def. 2.1.2 in [1]) 
%
% Syntax:  
%    val = supportFunc(E,dir)
%    val = supportFunc(E,dir,type)
%
% Inputs:
%    E   - ellipsoid object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the ellipsoid in the specified direction
%
% Example: 
%    E = ellipsoid([5 7;7 13],[1;2]);
%
%    val = supportFunc(E,[1;1]);
%   
%    figure
%    hold on
%    plot(E,[1,2],'b');
%    plot(conHyperplane(halfspace([1;1],val),[],[]),[1,2],'g');
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/supportFunc

% Author:       Victor Gassmann
% Written:      20-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'upper';
    
    if nargin >= 3 && ~isempty(varargin{1})
        type = varargin{1};
    end

    % compute support function
    l = dir/norm(dir);
    E0 = ellipsoid(E.Q);
    
    if strcmp(type,'upper')
        val = l'*center(E) + suppfnc(E0,l);
    else
        val = l'*center(E) - suppfnc(E0,l);
    end
    
    val = val * norm(dir);
    
end


% Auxiliary Functions -----------------------------------------------------

function [S] = suppfnc(E,l)

    S = zeros(size(l,2),1);
    
    for i=1:size(l,2)
        S(i) = l(:,i)'*E.q + sqrt(l(:,i)'*E.Q*l(:,i));
    end
end

%------------- END OF CODE --------------