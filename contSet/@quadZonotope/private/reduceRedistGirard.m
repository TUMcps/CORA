function qZ = reduceRedistGirard(qZ,order,options)
% reduceRedistGirard - Reduce remaining generators of a quadratic zonotope
% so that its order stays below a specified limit 
%
% Syntax:  
%    [qZred]=reduceRedistGirard(qZ,order,options)
%
% Inputs:
%    qZ - quadZonotope object
%    order - desired order of the zonotope
%    options - struct containing algorithm options
%
% Outputs:
%    qZ - reduced quadZonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      16-January-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if ~isempty(qZ.Grest)

        % determine dimension of zonotope
        dim=length(qZ.Grest(:,1));

        % only reduce if zonotope order is greater than the desired order
        if length(qZ.Grest(1,:))>dim*order

            % compute metric of generators (shortest generators)
            h=vnorm(qZ.Grest,1,2);
            [~,indices]=sort(h);

            % number of generators that are not reduced
            nUnreduced=max(0,floor(dim*(order-1)+1));
            
            % number of generators that are reduced
            nReduced=length(qZ.Grest(1,:))-nUnreduced;

            % pick generators that are reduced
            pickedGenerators=qZ.Grest(:,indices(1:nReduced));
            
            % scale generators in G for compensation
            [Gnew,rest] = reduction(qZ.G, pickedGenerators,options);
            qZ.G = Gnew;

            % unreduced generators
            qZ.Grest=[qZ.Grest(:,indices((nReduced+1):end)),rest];
        end
    end
end


function [E,rest] = reduction(E,G,options)
% substitute the generators in the matrix G by appropriate scaling of the
% generators in the matrix E

    % parse options struct
    if isfield(options,'reduceRedistGirard')
        if ~ismember(options.reduceRedistGirard,{'ascend','descend','none'})
           warning('Wrong value for "option.reduceRedistGirard"!. Only values "ascend", "descend" and "none" are supported.'); 
           opt = 'ascend';
        else
           opt = options.reduceRedistGirard;
        end
    else
        opt = 'ascend';
    end
    
    m = size(E,2);

    % compute the angles betwen the vectors in E and G
    ind = computeAngles(E,G);
    
    % assign independent generators to the closest dependent generator
    groups = cell(m,1);
    for i = 1:m
       groups{i} = find(ind == i); 
    end
    
    % construct new orthogonal basis for each dependent generator
    nBase = cell(m,1);
    
    for i = 1:m
        B = initialBase(E(:,i));
        nBase{i} = gramSchmidt(B);
    end
    
    % Heuristic: sort the generators by the projection error
    if ~strcmp(opt,'none')
        projErr = zeros(m,1);

        for i = 1:m
            % coordinate transformation to new space
            temp = nBase{i}' * G(:,groups{i});

            % compute projection error
            v = sum(abs(temp),2);
            projErr(i) = norm(v(2:end));      
        end

        % sort the dependent generators
        [~,index] = sort(projErr,opt);
        E = E(:,index);
        groups = groups(index);
        nBase = nBase(index);
    end
    
    % coordinate transformation and reduction with girard
    for i = 1:m
       if ~isempty(groups{i})
           % coordinate transformation
           N = nBase{i};
           G(:,groups{i}) = N' * G(:,groups{i});

           % reduce with girard
           v = sum(abs(G(:,groups{i})),2);
           E(:,i) = (1+abs(v(1))/norm(E(:,i)))*E(:,i);
           rest = diag(v);
           rest = N * rest(:,2:end);

           % assign the remaining generators to the remaining independent ones
           if i ~= m
               G = [G,rest];
               ind = computeAngles(E(:,i+1:end),rest)' + ones(1,size(rest,2))*i;

               for j = i+1:m
                  temp = find(ind == j);
                  temp = temp + ones(size(temp)) * (size(G,2)-size(rest,2));
                  groups{j} = [groups{j};temp'];
               end
           end
       end
    end
    
    % resort according to the original order
    if ~strcmp(opt,'none')
        E(:,index) = E;
    end
    
end


function B = initialBase(e)
% creates a initial basis that is later used to construct an orthonormal
% basis with the Gram-Schmidt method

    B = eye(length(e));
    e = e./norm(e);
    [~,ind] = max(e);
    B(:,ind) = B(:,1);
    B(:,1) = e;

end

function ind = computeAngles(E,G)
% finds for each independent generator (G) the index of the dependent
% generator (E) with the smallest angle 

    % compute the angles betwen the vectors
    temp = abs(G'*E);
    normE = sqrt(sum(E.^2,1));
    normG = sqrt(sum(G.^2,1));
    temp = temp .*(ones(size(temp,1),1)*(1./normE));
    temp = temp .*(ones(size(temp,2),1)*(1./normG))';

    % assign generators to the closest independent generator
    [~,ind] = max(temp,[],2);
end

function Q = gramSchmidt(B)
% constructs an orthonormal basis Q from the vectors in B

    Q = zeros(size(B));
    R = zeros(size(B));

    for j = 1:size(B,1)
        v = B(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)'*B(:,j);
            v = v-R(i,j)*Q(:,i);
        end
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
end