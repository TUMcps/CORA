function P = minus(P1,P2,varargin)
% minus - compute the Minkowski difference of two polytopes:
%         P1 - P2 = P <-> P + P2 \subseteq P1
%
% Syntax:  
%    P = minus(P1,P2)
%    P = minus(P1,P2,type)
%
% Inputs:
%    P1 - mptPolytope object
%    P2 - mptPolytope object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'approx')
%
% Outputs:
%    P - mptPolytope object after Minkowski difference
%
% Example: 
%    P1 = mptPolytope([1 0 -1 0 1;0 1 0 -1 1]',[4;4;4;4;4]);
%    P2 = mptPolytope([-1 0 1;-2 3 0]');
%
%    P = minus(P1,P2);
%
%    figure; hold on;
%    plot(P1);
%    plot(P2,[1,2],'r');
%    plot(P,[1,2],'g');
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
    if isnumeric(P2)
       P = P1 + (-P2); 
       
    else
        
        % parse input arguments
        type = 'exact';
        if nargin > 2 && ~isempty(varargin{1})
            type = varargin{1};
        end
        
        % exact computation 
        if strcmp(type,'exact')
        
            if isa(P2,'zonotope') && isa(P2,'interval')

                % convert to zonotope
                Z = zonotope(P2);

                % compute Minkowski diff. according to Theorem 1 in [1]
                c1 = center(Z);
                G = generators(Z);

                P = P1 + (-c1);

                for i = 1:size(G,2)
                    P = (P1 + G(:,i)) & (P1 + (-G(:,i)));
                end
                
            elseif isa(P2,'conZonotope') || isa(P2,'mptPolytope') || ...
                   isa(P2,'zonoBundle')

                % compute Minkowski difference according to Lemma 1 in [1]
                V = vertices(P2); P = P1;

                for i = 1:size(V,2)
                   P = P & (P1 + (-V(:,i))); 
                end
                
            else
                if strcmp(type,'exact')
                   error(errNoExactAlg(cZ1,cZ2));
                end
            end
        
        % under-approximative computation 
        elseif strcmp(type,'approx')
            
            % scale each halfspace of the polytope to obtain
            % inner-approximatin of the Minkowski difference
            A = P1.P.A;
            b = P1.P.b;
            
            for i = 1:size(A,1)
               l = supportFunc(P2,A(i,:)','upper');
               b(i) = b(i) - l;
            end
            
            P = mptPolytope(A,b);          
            
        else
            [id,msg] = errWrongInput('type');
            error(id,msg);       
        end
    end
end

%------------- END OF CODE --------------