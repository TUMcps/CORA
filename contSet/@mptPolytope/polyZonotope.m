function pZ = polyZonotope(obj)
% polyZonotope - Convert polytope to a polynomial zonotope
%
% Syntax:  
%    pZ = polyZonotope(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    A = [0 -1;2 -1;1 2;-5 1;-1 -1];
%    b = [2;4;7;9;3];
%    poly = mptPolytope(A,b);
%
%    pZ = polyZonotope(poly);
%
%    plot(poly,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-3,4]);
%    ylim([-3,5]);
% 
%    figure
%    plot(pZ,[1,2],'b','Filled',true,'EdgeColor','none','Splits',12);
%    xlim([-3,4]);
%    ylim([-3,5]);   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, taylm/polyZonotope

% Author:       Niklas Kochdumper
% Written:      26-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% distinguish between 2D and multi-dimensional case
if size(obj.P.A,2) == 2
 
    % compute polytope vertices
    V = vertices(obj);
    
    if size(V,2) > 2
        ind = convhull(V(1,:),V(2,:));
        V = V(:,ind);
    end
    
    % loop over all polytope faces
    pZface = cell(ceil(size(V,2)/2),1);
    counter = 1;
    i = 1;
    
    while i < size(V,2)
        
        % construct polynomial zonotope from the facet
        c = 0.5*(V(:,i) + V(:,i+1));
        g = 0.5*(V(:,i) - V(:,i+1));
        
        pZface{counter} = polyZonotope(c,g,[],1);
        
        counter = counter + 1;
        i = i + 2;
    end
    
    % consider the last face
    if isempty(pZface{end})
        c = 0.5*(V(:,end-1) + V(:,end));
        g = 0.5*(V(:,end-1) - V(:,end));      
        pZface{counter} = polyZonotope(c,g,[],1);
    end
    
    % iteratively compute the convex hull of neightbouring faces until only
    % one single polynomial zonotope is left
    while length(pZface) > 1
       
        pZtemp = cell(floor(length(pZface)/2),1);
        
        counter = 1;
        counterNew = 1;
        
        % loop over all pairs
        while counter+1 <= length(pZface)
           
            % construct polyZonotope objects with appropriate id vectors
            pZ1_ = pZface{counter}; 
            pZ2_ = pZface{counter+1};
            
            id1 = 1:length(pZ1_.id);
            id2 = length(pZ1_.id)+1 : length(pZ1_.id) + length(pZ2_.id);
            id2(end) = id1(end);
            
            G2 = pZ2_.G;
            ind = find(pZ2_.expMat(end,:) > 0);
            G2(:,ind) = -G2(:,ind);

            pZ1 = polyZonotope(pZ1_.c,pZ1_.G,pZ1_.Grest,pZ1_.expMat,id1');
            pZ2 = polyZonotope(pZ2_.c,G2,pZ2_.Grest,pZ2_.expMat,id2'); 
            
            % compute convex hull
            pZtemp{counterNew} = enclose(pZ1,pZ2);
            counterNew = counterNew + 1;
            counter = counter + 2;     
        end
        
        if counterNew == ceil(length(pZface)/2)
            pZtemp = [pZtemp; pZface(end)]; 
        end
        
        pZface = pZtemp;      
    end
    
    pZ = pZface{1};
    
else
    
    % compute polytope vertices
    V = vertices(obj);
    
    % convert each vertex to a polynomial zonotope
    list = cell(size(V,2),1);
    cent = V;
    lenID = zeros(length(list),1);
    
    for i = 1:length(list)
       list{i} = polyZonotope(V(:,i),[],[],[]); 
    end
    
    % recursively compute the convex hull
    while length(list) > 1
        
        tempList = cell(length(list),1);
        tempCent = zeros(size(cent));
        tempLenID = zeros(length(list),1);
        indID = find(lenID == max(lenID));
        counter = 1;
        
        % loop over all polynomial zonotopes in the list
        while ~isempty(list) && length(list) > 1
            
            % determine nearest polynomial zonotope from the list
            if length(indID) >= 2
                ind = nearestNeighbour(cent(:,indID(1)),cent(:,indID(2:end)));
                ind1 = indID(1);
                ind2 = indID(ind+1);
                indID = [indID(2:ind)-1;indID(ind+2:end)-2];
            else
                ind1 = 1;
                ind2 = nearestNeighbour(cent(:,1),cent(:,2:end)) + 1;
            end
            
            % construct polyZonotope objects with appropriate id vectors
            pZ1_ = list{ind1}; 
            pZ2_ = list{ind2};
            
            if ~isempty(pZ1_.expMat) && ~isempty(pZ2_.expMat) 
                id1 = 1:length(pZ1_.id);
                id2 = length(pZ1_.id)+1 : length(pZ1_.id) + length(pZ2_.id);

                pZ1 = polyZonotope(pZ1_.c,pZ1_.G,pZ1_.Grest,pZ1_.expMat,id1');
                pZ2 = polyZonotope(pZ2_.c,pZ2_.G,pZ2_.Grest,pZ2_.expMat,id2');            
            else
                pZ1 = pZ1_;
                pZ2 = pZ2_;
            end
            
            % compute convex hull
            tempList{counter} = enclose(pZ1,pZ2);
            tempCent(:,counter) = center(tempList{counter});
            
            % update id-vectors if empty
            if isempty(tempList{counter}.id)
               temp = tempList{counter};
               tempList{counter} = polyZonotope(temp.c,temp.G, ...
                                                temp.Grest,temp.expMat,1);
            end
            
            % update length of ID vectors
            tempLenID(counter) = length(tempList{counter}.id);
            
            % update variables
            cent(:,[ind1,ind2]) = [];
            list([ind1,ind2]) = [];
            counter = counter + 1;
        end
        
        % add last list element 
        if ~isempty(list)
            tempList{counter} = list{1};
            tempCent(:,counter) = cent(:,1);
            tempLenID(counter) = length(list{1}.id);
        else
            counter = counter - 1;
        end
        
        % update lists
        list = tempList(1:counter);
        cent = tempCent(:,1:counter); 
        lenID = tempLenID(1:counter);
    end

    % construct the resulting polynomial zonotopes
    pZ = list{1};
end 
    
    
% Auxiliary Functions -----------------------------------------------------

function ind = nearestNeighbour(p,points)

    dist = sum((points-p*ones(1,size(points,2))).^2,1);
    [~,ind] = min(dist);
    ind = ind(1);

%------------- END OF CODE --------------