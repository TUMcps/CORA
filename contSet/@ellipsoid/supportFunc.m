function [val,x] = supportFunc(E,dir,varargin)
% supportFunc - Calculate the upper or lower bound of an ellipsoid object
%               along a certain direction (see Def. 2.1.2 in [1]) 
%
% Syntax:  
%    [val,x] = supportFunc(E,dir)
%    [val,x] = supportFunc(E,dir,type)
%
% Inputs:
%    E   - ellipsoid object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the ellipsoid in the specified direction
%    x   - point for which holds: dir'*x=val
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
% Last update:  12-March-2021
%               27-July-2021 (fixed degenerate case)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'upper';
    
    if nargin >= 3 && ~isempty(varargin{1})
        type = varargin{1};
    end

    if strcmp(type,'upper')
        s = 1;
    else
        s = -1;
        
    end
    val = dir'*E.q + s*sqrt(dir'*E.Q*dir);
    if nargout > 1
        H = conHyperplane(dir',val);
        n_nd = rank(E);
        n = dim(E);
        T = eye(n);
        qt_rem = [];
        
        if rank(E)==0
            x = E.q;
            return;
        elseif n_nd<n
            % remove degenerate dimension and reduce E and hyperplane
            [T,~,~] = svd(E.Q);
            dir_t = dir'*T;
            
            % check if chosen direction is aligned with one of the
            % degenerate directions
            if any(withinTol(dir_t(1:n_nd),0,E.TOL))
                % this results in a 0=0 type equation (so no constraint on
                % the remaining variables) => choose center as point
                x = E.q;
                return;
            end
            E = T'*E;
            qt_rem = E.q(n_nd+1:end);
            % project both into non-degenerate dimensions
            H = conHyperplane(dir_t(1:n_nd),val-dir_t(n_nd+1:end)*qt_rem);
            E = project(E,1:n_nd);
        end
        % find point x which attains the value by forming hyperplane and
        % intersecting with ellipsoid
        E_res = E&H;
        assert(~isempty(E_res),'intersection with tangent hp should not be empty');
        % basically contains only 1 point
        x = T*[E_res.q;qt_rem];
    end
end
%------------- END OF CODE --------------