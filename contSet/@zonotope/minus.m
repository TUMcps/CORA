function [Z] = minus(minuend,subtrahend)
% minus - Overloaded '-' operator for approximating the Minkowski difference of two
% zonotopes or a zonotope with a vector. A - B = C <-> B + C \subseteq A
%
% Syntax:  
%    [Z] = minus(minuend,subtrahend)
%
% Inputs:
%    minuend - zonotope object
%    subtrahend - zonotope object or numerical vector
%
% Outputs:
%    Z - Zonotpe after Minkowsi difference
%
% Example: 
%    Z1 = zonotope([1 2 2; 0 0 2]);
%    Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
%    Z3 = Z1 - Z2;
%    Z4 = Z2 + Z3;
%    plot(Z1,[1 2], 'b');
%    hold on
%    plot(Z2,[1 2], 'r');
%    plot(Z3,[1 2], 'g');
%    plot(Z4,[1 2], 'k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      10-June-2015
% Last update:  22-July-2015
%               05-August-2015
%               20-August-2015
%               30-July-2016
%               14-November-2016
% Last revision:---

%------------- BEGIN CODE --------------


%Is subtrahend a zonotope?
if isa(subtrahend,'zonotope')

    
    %determine generators to be kept---------------------------------------
    %obtain halfspace representation
    [P,comb] = polytope(minuend);
    
%     %test how many halfspaces are removed before applying Minkowski
%     %difference
%     %remove redundant halfspaces
%     Pnew = reduce(P);
    
    P_h = halfspace(P);
    HorigTwice = get(P_h,'H');
    KorigTwice = get(P_h,'K');
    Horig = HorigTwice(1:0.5*end,:);
    
    %nr of subtrahend generators
    subtrahendGens = length(subtrahend.Z(1,:)) - 1;
    
    %intersect polytopes
    delta_K = 0*KorigTwice;
    for i = 1:subtrahendGens
        delta_K = delta_K + abs(HorigTwice*subtrahend.Z(:,i+1));
    end
    Korig_new = KorigTwice - delta_K;
    P_int = mptPolytope(HorigTwice,Korig_new);
    %P_int = mptPolytope(Horig,Korig_new(1:0.5*end));
    
    %remove redundant halfspaces and remember inices
    removedInd = removedHalfspaces(P_int,Horig,Korig_new(1:0.5*end));
    %removedIndAlt = removedHalfspacesAlternative(P_int,Horig,Korig_new(1:0.5*end));
    
    if removedInd == -inf
        % Minkowski difference does not exist
        Z = [];
    else
        if ~isempty(removedInd)
            %count generators that have been removed
            gens = length(minuend.Z(1,2:end));
            indices = zeros(1,gens);
            for i=1:length(removedInd)
                addIndices = zeros(1,gens);
                addIndices(comb(removedInd(i),:)) = 1;
                indices = indices + addIndices;
            end
            % find generators to be removed
            n = length(minuend.Z(:,1));
            requiredRemovals = factorial(gens-1)/(factorial(n-2)*factorial(gens-n+1));  %binom(gens-1,dim-2);
            indRemove = (indices == requiredRemovals);

            %obtain reduced minuend
            G = minuend.Z(:,2:end);
            G(:,indRemove) = [];
            minuend.Z = [minuend.Z(:,1), G];

            %remove H, K
            C = Horig;
            C(removedInd,:) = [];
            d = Korig_new(1:0.5*end,:);
            d(removedInd) = [];
        else
            C = Horig;
            d = Korig_new(1:0.5*end,:);
        end
    
        %compute center
        c = minuend.Z(:,1) - subtrahend.Z(:,1);

        %obtain minuend generators
        G = minuend.Z(:,2:end);
        %----------------------------------------------------------------------

        %reverse computation from halfspace generation
        delta_d = d - C*minuend.Z(:,1);
        A_abs = abs(C*G);
        %alpha = A_abs\delta_d; %solve linear set of equations
        alpha = pinv(A_abs)*delta_d; %solve linear set of equations
        for i=1:length(alpha)
            Gnew(:,i) = alpha(i)*minuend.Z(:,i+1);
        end

        % instantiate Z
        Z = zonotope([c,Gnew]); 
    end

    
%is subtrahend a vector?
elseif isnumeric(subtrahend)
    %instantiate Z
    Z = minuend;
    
    %Calculate minkowski difference
    Z.Z(:,1)=Z.Z(:,1) - subtrahend;
end

% %Is subtrahend a zonotope?
% if isa(subtrahend,'zonotope')
% 
%     
%     %determine generators to be kept---------------------------------------
%     %obtain halfspace representation
%     [P,comb] = polytope(minuend);
%     P_h = halfspace(P);
%     HorigTwice = get(P_h,'H');
%     KorigTwice = get(P_h,'K');
%     Horig = HorigTwice(1:0.5*end,:);
%     
%     %nr of subtrahend generators
%     subtrahendGens = length(subtrahend.Z(1,:)) - 1;
%     
%     %intersect polytopes
%     delta_K = 0*KorigTwice;
%     for i = 1:subtrahendGens
%         delta_K = delta_K + abs(HorigTwice*subtrahend.Z(:,i+1));
%     end
%     Korig_new = KorigTwice - delta_K;
%     P_int = mptPolytope(HorigTwice,Korig_new);
%     
%     if isempty(P_int)
%         % Minkowski difference does not exist
%         Z = [];
%     else
%         
%         %remove H, K
%         C = Horig;
%         d = Korig_new(1:0.5*end,:);
%     
%         %compute center
%         c = minuend.Z(:,1) - subtrahend.Z(:,1);
% 
%         %obtain minuend generators
%         G = minuend.Z(:,2:end);
%         %----------------------------------------------------------------------
% 
%         %reverse computation from halfspace generation
%         delta_d = d - C*minuend.Z(:,1);
%         A_abs = abs(C*G);
%         alpha = A_abs\delta_d; %solve linear set of equations
%         for i=1:length(alpha)
%             Gnew(:,i) = alpha(i)*minuend.Z(:,i+1);
%         end
% 
%         % instantiate Z
%         Z = zonotope([c,Gnew]); 
%     end
% 
%     
% %is subtrahend a vector?
% elseif isnumeric(subtrahend)
%     %instantiate Z
%     Z = minuend;
%     
%     %Calculate minkowski difference
%     Z.Z(:,1)=Z.Z(:,1) - subtrahend;
% end


%------------- END OF CODE --------------