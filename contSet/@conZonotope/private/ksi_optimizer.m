function Dksi = ksi_optimizer(obj)    
% ksi_optimizer - determine the tighend domains for the zonotope factors
%                 ksi by solving a linear program
%
% Syntax:  
%    Dksi = ksi_optimizer(obj)
%
% Inputs:
%    obj - constrained zonotope object
%
% Outputs:
%    Dksi - new tighend domains for the zonotope factors (class: interval)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if ~isempty(obj.A)
        
        % object properties  
        n = size(obj.A, 2);
        A = obj.A;
        b = obj.b;

        % ksi in [-1, 1]
        lb = -ones(n,1);
        ub = ones(n,1);

        % initialize ksi
        ksi_min = zeros(n,n);
        ksi_max = zeros(n,n);

        % options
        options = optimoptions('linprog','Algorithm','dual-simplex', 'display','off');

        for i = 1:n

            % min/max ksi_i
            f = zeros(n, 1);
            f(i, 1) = 1;

            % minimize (Equation (25) in [1])
            [x, ~, flag_min] = linprog(f,[],[],A,b,lb,ub,options);
            ksi_min(:,i) = x;

            % maximize  (Equation (26) in [1])
            [x, ~, flag_max] = linprog(-f,[],[],A,b,lb,ub,options);
            ksi_max(:,i) = x;
            if flag_min ~= 1 || flag_max ~= 1
                error('Optimization error')
            end
        end

        % delete dublicates
        ksi = unique([ksi_min, ksi_max]','rows');
        ksi = ksi';

        % calculate tightend domain for the zonotope factors
        ksi_l = min(ksi,[],2);
        ksi_u = max(ksi,[],2);

        Dksi = interval(ksi_l,ksi_u);


    else        % no constraints -> return unit-cube as domain for ksi

        n = size(obj.Z,2)-1;
        Dksi = interval(-ones(n,1),ones(n,1));
    end
    
end   

%------------- END OF CODE --------------