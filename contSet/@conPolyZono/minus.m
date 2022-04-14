function cPZ = minus(cPZ1,cPZ2,varargin)
% minus - compute the Minkowski difference of two conPolyZono objects:
%         cPZ1 - cPZ2 = cPZ <-> cPZ + cPZ2 \subseteq cPZ1
%
% Syntax:  
%    cPZ = minus(cPZ1,cPZ2)
%    cPZ = minus(cPZ1,cPZ2,type)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'approx')
%
% Outputs:
%    P - conPolyZono object after Minkowski difference
%
% Example: 
%   cPZ = conPolyZono([3;3],[1 -2 1; 2 3 1],[1 0 2;0 1 1]);
%   zono = zonotope([0;0],[0.4;0.4]);
%
%   res = cPZ - zono;
%
%   figure; hold on;
%   plot(cPZ,[1,2],'b');
%   plot(res,[1,2],'r','Splits',25);
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minus, conZonotope/minus

% Author:       Niklas Kochdumper
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % different algorithms for different set representations
    if isnumeric(cPZ2)
       cPZ = cPZ1 + (-cPZ2); 
       
    else
        
        % consider independent generators
        if nargin > 2 && strcmp(varargin{1},'exact') && ~isempty(cPZ1.Grest)
            cPZ1 = removeIndepGens(cPZ1);
        else
            cPZ1.Grest = [];
        end
        
        % different cases for different set representations
        if isa(cPZ2,'zonotope') || isa(cPZ2,'interval')

            % convert to zonotope
            Z = zonotope(cPZ2);

            % compute Minkowski difference according to Theorem 1 in [1]
            c = center(Z);
            G = generators(Z);

            cPZ = cPZ1 + (-c);

            for i = 1:size(G,2)
                cPZ = (cPZ1 + G(:,i)) & (cPZ1 + (-G(:,i)));
            end

        elseif isa(cPZ2,'conZonotope') || isa(cPZ2,'mptPolytope') || ...
               isa(cPZ2,'zonoBundle')

           % parse input arguments
            type = 'exact';
            if nargin > 2 && ~isempty(varargin{1})
                type = varargin{1};
            end
           
            % compute exact result or fast inner-approximation
            if strcmp(type,'approx')
               cPZ = minus(cPZ2,zonotope(cPZ2));
               return;
            end

            % compute exact result according to Lemma 1 in [1]
            V = vertices(cPZ2); cPZ = cPZ1;

            for i = 1:size(V,2)
               cPZ = cPZ & (cPZ1 + (-V(:,i))); 
            end

        else

            % parse input arguments
            type = 'approx';
            if nargin > 2 && ~isempty(varargin{1})
                type = varargin{1};
            end
            
            if strcmp(type,'exact')
               error(errNoExactAlg(cPZ1,cPZ2));
            end

            % compute inner-approximation by enclosing second set with zonotope
            cPZ = cPZ1 - zonotope(cPZ2);
        end
    end 
end


% Auxiliary Functions -----------------------------------------------------

function cPZ = removeIndepGens(cPZ)
% redefine independent generators as new dependent generators

    E = eye(size(cPZ.Grest,2));
    id = [cPZ.id;(max(cPZ.id)+1:max(cPZ.id+size(cPZ.Grest,2)))'];
    if isempty(cPZ.A)
       cPZ = conPolyZono(cPZ.c,[cPZ.G,cPZ.Grest], ...
                          blkdiag(cPZ.expMat,E),[],id);
    else
       cPZ = conPolyZono(cPZ.c,[cPZ.G,cPZ.Grest], ...
                          blkdiag(cPZ.expMat,E),cPZ.A,cPZ.b, ...
                          blkdiag(cPZ.expMat_,E),[],id);
    end
end

%------------- END OF CODE --------------