function res = isPolytope(pZ)
% isPolytope - Checks if a polynomial zonotope represents a polytope
%
% Syntax:
%    res = isPolytope(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([-0.5;0],[1.5 -0.5 -0.5;-0.5 -2 0.5],[],[1 0 1;0 1 1]);
%    pZ2 = polyZonotope([-0.5;0],[-0.5 -0.5 1.5;0.5 -2 -0.5],[],[1 0 1;0 1 1]);
%   
%    isPolytope(pZ1)
%    isPolytope(pZ2)
%
%    figure; hold on
%    plot(pZ1,[1,2],'FaceColor','b','Splits',10);
%
%    figure; hold on
%    plot(pZ2,[1,2],'FaceColor','r','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope

% Authors:       Niklas Kochdumper
% Written:       08-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = true;
    
    % check if variable appears with exponent greater than 1
    if any(any(pZ.expMat >= 2))
       res = false;
       return;
    end

    % compute vertices
    [V,I] = aux_polyVertices(pZ);

    % loop over all vertices
    for i = 1:size(V,2)
       
        % compute vectors of normal cone
        N = aux_normalCone(pZ,I{i});
        
        % loop over all vertices
        for j = 1:size(V,2)
           if i ~= j
               if ~aux_inCone(N,V(:,i),V(:,j))
                   res = false;
                   return;
               end
           end
        end      
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_inCone(N,v,v_)
% checks if a vertex is inside the normal cone

    % construct constraints and objective function
    [n,m] = size(N);

    A1 = [N,-eye(n)];
    b1 = v_-v;
    A2 = [-N,-eye(n)];
    b2 = v-v_;
    
    lb = zeros(n+m,1);
    
    f = [zeros(1,m),ones(1,n)];
    
    % solve linear program
    persistent options
    if isempty(options)
        options = optimoptions('linprog','Display','off');
    end
    
    x = linprog(f,[A1;A2],[b1;b2],[],[],lb,[],options);
    
    % check if problem is solvable
    res = all(x(m+1:end) < eps);

end

function Nall = aux_normalCone(pZ,Ilist)
% compute the normal cone for a vertex
    
    % loop over all factor combinations resulting in the same vertex
    Nall = [];
    
    for l = 1:size(Ilist,2)
        
        I = Ilist(:,l);
        N = zeros(size(pZ.G,1),length(I));

        % loop over all factors
        for i = 1:length(I)

            ind = find(pZ.expMat(i,:) > 0);

            % loop over all generators
            for j = 1:length(ind)
               k = ind(j);
               N(:,i) = N(:,i) + pZ.G(:,k) * prod(I.^pZ.expMat(:,k))/I(i)^pZ.expMat(i,k);
            end
        end
    
        N = N.*(-sign(I)');
        N = N(:,sum(abs(N),1) > eps);
        Nall = [Nall,N];
    end
    
end

function [V,Ilist] = aux_polyVertices(pZ)
% compute the polytope vertices and store the corresponding factor values

    % determine all potential vertices
    p = size(pZ.expMat,1);
    n = size(pZ.G,1);
    
    I = vertices(interval(-ones(p,1),ones(p,1)));
    
    V = zeros(n,size(I,2));
    
    for i = 1:size(I,2)
       V(:,i) = pZ.c;
       for j = 1:size(pZ.G,2)
           V(:,i) = V(:,i) + pZ.G(:,j) * prod(I(:,i).^pZ.expMat(:,j));
       end
    end
    
    % remove redundant points
    [V,ind] = sortrows(V');
    
    V = V';
    I = I(:,ind);

    V_ = zeros(size(V));
    Ilist = cell(size(V,2),1);
    
    V_(:,1) = V(:,1);
    Itemp = I(:,1);
    counter = 1;

    for i = 2:size(V,2)
        if ~all(abs(V(:,i)-V_(:,counter)) < 1e-10)
           counter = counter + 1;
           V_(:,counter) = V(:,i);
           Ilist{counter-1} = Itemp;
           Itemp = I(:,i);
        else
           Itemp = [Itemp,I(:,i)];
        end
    end
    
    Ilist{counter} = Itemp;

    V = V_(:,1:counter);
    Ilist = Ilist(1:counter);

    % determine vertices with the n-dimensional convex hull
    ind = convhulln(V');
    ind = unique(squeeze(ind));
    
    V = V(:,ind);
    Ilist = Ilist(ind);

end

% ------------------------------ END OF CODE ------------------------------
