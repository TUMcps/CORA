function cZ = minus(cZ1,cZ2,varargin)
% minus - compute the Minkowski difference of two constrained zonotopes:
%         cZ1 - cZ2 = cZ <-> cZ + cZ2 \subseteq cZ1
%
% Syntax:  
%    cZ = minus(cZ1,cZ2)
%    cZ = minus(cZ1,cZ2,type)
%
% Inputs:
%    cZ1 - conZonotope object
%    cZ2 - conZonotope object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'approx')
%
% Outputs:
%    cZ - conZonotope object after Minkowski difference
%
% Example: 
%    Z = [0 4 0 0;0 0 4 0];
%    A = [1 1 1.5];
%    b = -0.5;
%    cZ1 = conZonotope(Z,A,b);
%
%    Z = [1 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZ2 = conZonotope(Z,A,b);
%
%    cZ = minus(cZ1,cZ2);
%
%    figure; hold on;
%    plot(cZ1);
%    plot(cZ2,[1,2],'r');
%    plot(cZ,[1,2],'g');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minus

% Author:       Niklas Kochdumper
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % different algorithms for different set representations
    if isnumeric(cZ2)
       cZ = cZ1 + (-cZ2); 
       
    elseif isa(cZ2,'zonotope') && isa(cZ2,'interval')
        
        % convert to zonotope
        Z = zonotope(cZ2);
        
        % compute Minkowski difference according to Theorem 1 in [1]
        c = center(Z);
        G = generators(Z);

        cZ = cZ1 + (-c);

        for i = 1:size(G,2)
            cZ = (cZ1 + G(:,i)) & (cZ1 + (-G(:,i)));
        end
        
    elseif isa(cZ2,'conZonotope') || isa(cZ2,'mptPolytope') || ...
           isa(cZ2,'zonoBundle')
       
        % parse input arguments
        type = 'exact';
        if nargin > 2 && ~isempty(varargin{1})
            type = varargin{1};
        end
        
        % compute exact result or fast inner-approximation
        if strcmp(type,'approx')
           cZ = minus(cZ2,zonotope(cZ2));
           return;
        end
        
        % compute exact result according to Lemma 1 in [1]
        V = vertices(cZ2); cZ = cZ1;
        
        for i = 1:size(V,2)
           cZ = cZ & (cZ1 - V(:,i)); 
        end
        
    else
        
        % parse input arguments
        type = 'exact';
        if nargin > 2 && ~isempty(varargin{1})
            type = varargin{1};
        end
        
        if strcmp(type,'exact')
           error(errNoExactAlg(cZ1,cZ2));
        end
        
        % compute inner-approximation by enclosing second set with zonotope
        cZ = cZ1 - zonotope(cZ2);
    end
end

%------------- END OF CODE --------------