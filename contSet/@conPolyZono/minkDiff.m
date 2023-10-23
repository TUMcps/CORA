function cPZ = minkDiff(cPZ1,cPZ2,varargin)
% minkDiff - compute the Minkowski difference of two conPolyZono objects:
%         cPZ1 - cPZ2 = cPZ <-> cPZ + cPZ2 \subseteq cPZ1
%
% Syntax:
%    cPZ = minkDiff(cPZ1,cPZ2)
%    cPZ = minkDiff(cPZ1,cPZ2,type)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'inner')
%
% Outputs:
%    P - conPolyZono object
%
% Example: 
%   cPZ = conPolyZono([3;3],[1 -2 1; 2 3 1],[1 0 2;0 1 1]);
%   Z = zonotope([0;0],[0.4;0.4]);
%
%   res = minkDiff(cPZ,Z);
%
%   figure; hold on;
%   plot(cPZ,[1,2],'b');
%   plot(res,[1,2],'r','Splits',12);
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff, conZonotope/minkDiff

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   09-November-2022 (MW, rename 'minkDiff')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    type = setDefaultValues({'exact'},varargin);

    % check input arguments
    inputArgsCheck({{cPZ1,'att','conPolyZono'};
                    {cPZ2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','inner'}}});
    
    % different algorithms for different set representations
    if isnumeric(cPZ2)
       cPZ = cPZ1 + (-cPZ2); 
       
    else
        
        % consider independent generators
        if strcmp(type,'exact') && ~isempty(cPZ1.GI)
            cPZ1 = aux_removeIndepGens(cPZ1);
        else
            cPZ1.GI = [];
        end
        
        % different cases for different set representations
        if isa(cPZ2,'zonotope') || isa(cPZ2,'interval')

            % convert to zonotope
            Z = zonotope(cPZ2);

            % compute Minkowski difference according to Theorem 1 in [1]
            c = Z.c;
            G = Z.G;

            cPZ = cPZ1 + (-c);

            for i = 1:size(G,2)
                cPZ = and_(cPZ + G(:,i),cPZ + (-G(:,i)),'exact');
            end

        elseif isa(cPZ2,'conZonotope') || isa(cPZ2,'polytope') || ...
               isa(cPZ2,'zonoBundle')
           
            % compute exact result or fast inner-approximation
            if strcmp(type,'inner')
                cPZ = minkDiff(cPZ2,zonotope(cPZ2));
                return;
            end

            % compute exact result according to Lemma 1 in [1]
            V = vertices(cPZ2);
            cPZ = cPZ1 + (-V(:,1));

            for i = 2:size(V,2)
                cPZ = and_(cPZ,cPZ1 + (-V(:,i)),'exact'); 
            end

        else
            
            if strcmp(type,'exact')
                throw(CORAerror('CORA:noExactAlg',cPZ1,cPZ2));
            end

            % compute inner-approximation by enclosing second set with zonotope
            cPZ = minkDiff(cPZ1,zonotope(cPZ2));
        end
    end

end


% Auxiliary functions -----------------------------------------------------

function cPZ = aux_removeIndepGens(cPZ)
% redefine independent generators as new dependent generators

    E = eye(size(cPZ.GI,2));
    id = [cPZ.id;(max(cPZ.id)+1:max(cPZ.id+size(cPZ.GI,2)))'];
    if isempty(cPZ.A)
       cPZ = conPolyZono(cPZ.c,[cPZ.G,cPZ.GI], blkdiag(cPZ.E,E),[],id);
    else
       cPZ = conPolyZono(cPZ.c,[cPZ.G,cPZ.GI], blkdiag(cPZ.E,E), ...
           [cPZ.A,zeros(1,size(cPZ.GI,2))],cPZ.b, blkdiag(cPZ.EC,E),[],id);
    end
end

% ------------------------------ END OF CODE ------------------------------
