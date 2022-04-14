function res = cubMap(cZ,varargin)
% cubMap - computes an enclosure of the set corresponding to the cubic 
%          multiplication of a constrained zonotope with a third-order tensor
%
% Description:
%   Calculates the following set:
%   { z = (x' T x) * x | x \in cZ }
%
%   If three conZonotopes are provided, the function calculates the set:
%   { z = (x1' T x2) * x3 | x1 \in cZ1, x2 \in cZ2, x3 \in cZ3 }
%
% Syntax:  
%    res = cubMap(cZ,T)
%    res = cubMap(cZ1,cZ2,cZ3,T)
%    res = cubMap(cZ,T,ind)     
%    res = cubMap(cZ1,cZ2,cZ3,T,ind)
%
% Inputs:
%    cZ,cZ1,cZ2,cZ3 - conZonotope objects
%    T - third-order tensor
%    ind - cell-array containing the non-zero indizes of the tensor
%
% Outputs:
%    res - conZonotope object representing the set of the cubic mapping
%
% Example: 
%    % cubic multiplication
%    cZ = conZonotope.generateRandom(2,[],4);
%    
%    T{1,1} = rand(2);
%    T{1,2} = rand(2);
%    T{2,1} = rand(2);
%    T{2,2} = rand(2);
%
%    cZcub = cubMap(cZ,T);
%
%    figure 
%    subplot(1,2,1)
%    plot(cZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    subplot(1,2,2)
%    plot(cZcub,[1,2],'b','Filled',true,'EdgeColor','none','Template',50);
%
%    % mixed cubic multiplication
%    Z2 = zonotope.generateRandom(2,[],5);
%    Z3 = zonotope.generateRandom(2,[],3);
%
%    ZcubMixed = cubMap(Z,Z2,Z3,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadMap, zonotope/cubMap

% Author:       Niklas Kochdumper
% Written:      02-November-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % cubic multiplication or mixed cubic multiplication
    if isa(varargin{1},'conZonotope')
        
        % check user input
        if nargin < 4
           error('Wrong syntax for function cubMap!'); 
        end
        
        cZ2 = varargin{1};
        cZ3 = varargin{2};
        T = varargin{3};
        
        % parse optional input arguments
        if nargin > 4
           ind = varargin{4}; 
        else
           temp = 1:size(T,2);
           ind = repmat({temp},[size(T,1),1]);
        end 
        
        % mixed cubic multiplication
        res = cubMapMixed(cZ,cZ2,cZ3,T,ind);
        
    else
        
        % check user input
        if nargin < 2
           error('Wrong syntax for function cubMap!'); 
        end
        
        T = varargin{1};
        
        % parse optional input arguments
        if nargin > 2
           ind = varargin{2}; 
        else
           temp = 1:size(T,2);
           ind = repmat({temp},[size(T,1),1]);
        end 
        
        % cubic multiplication
        res = cubMapSingle(cZ,T,ind);       
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = cubMapSingle(cZ,T,ind)
% calulates the following set:      { z = (x' T x) * x | x \in cZ }

    if isempty(cZ.A) 

        temp = cubMap(zonotope(cZ.Z),T,ind);
        res = conZonotope(temp);
        
    else
        
        % rescale constrained zonotope to reduce over-approximation
        if ~isempty(cZ.A)
            cZ = rescale(cZ);
        end
        
        % compute cubic map for polynomial zonotope to keep track of
        % dependencies
        pZ = polyZonotope(zonotope(cZ.Z));
        temp = cubMap(pZ,T,ind);

        % compute constraint matrix from con. mat. of original conZonotope
        index = find(sum(temp.expMat,1) == 1);
        A = zeros(size(cZ.A,1),size(temp.G,2));
        b = cZ.b;

        for i = 1:size(temp.expMat,1)
            ind_ = find(temp.id == i);
            if ~isempty(ind_)
                ind = find(temp.expMat(ind_,index) == 1);
                if ~isempty(ind)
                   A(:,index(ind)) = cZ.A(:,i); 
                end
            end
        end

        % compute zonotope enclosure of the polynomial zonotopes
        ind = find(prod(ones(size(temp.expMat))-mod(temp.expMat,2),1) == 1);

        c = temp.c + 0.5 * sum(temp.G(:,ind),2);
        G = temp.G;
        G(:,ind) = 0.5*G(:,ind);
        
        % remove all-zero constraints
        ind = find(sum(abs(A),2) == 0);
        A(ind,:) = [];
        b(ind) = [];

        % construct resulting constrained zonotope
        res = conZonotope(c,G,A,b);
    end
    
end

function res = cubMapMixed(cZ1,cZ2,cZ3,T,ind)
% calulates the following set:
% { z = (x1' T x2) * x3 | x1 \in cZ1, x2 \in cZ2, x3 \in cZ3 }

    % rescale constrained zonotope to reduce over-approximation
    if ~isempty(cZ1.A)
        cZ1 = rescale(cZ1);
    end
    if ~isempty(cZ2.A)
        cZ2 = rescale(cZ2);
    end
    if ~isempty(cZ3.A)
        cZ3 = rescale(cZ3);
    end

    % compute cubic map for polynomial zonotope to keep track of
    % dependencies
    m1 = size(cZ1.Z,2) - 1; ind1 = 1:m1;
    m2 = size(cZ2.Z,2) - 1; ind2 = m1+1:m1+m2;
    m3 = size(cZ3.Z,2) - 1; ind3 = m1+m2+1:m1+m2+m3;
    
    pZ1 = polyZonotope(cZ1.Z(:,1),cZ1.Z(:,2:end),[],eye(m1),ind1');
    pZ2 = polyZonotope(cZ2.Z(:,1),cZ2.Z(:,2:end),[],eye(m2),ind2');
    pZ3 = polyZonotope(cZ3.Z(:,1),cZ3.Z(:,2:end),[],eye(m3),ind3');
    
    temp = cubMap(pZ1,pZ2,pZ3,T,ind);

    % compute constraint matrix from con. mat. of original conZonotope
    index = find(sum(temp.expMat,1) == 1);
    A = []; b = [];
    
    list = {cZ1,cZ2,cZ3};
    listInd = {ind1,ind2,ind3};
    
    for j = 1:length(list)
    
        cZ = list{j};
        
        if ~isempty(cZ.A)
            
            Atemp = zeros(size(cZ.A,1),size(temp.G,2));

            for i = listInd{j}
                ind_ = find(temp.id == i);
                if ~isempty(ind_)
                    ind = find(temp.expMat(ind_,index) == 1);
                    if ~isempty(ind)
                       Atemp(:,index(ind)) = cZ.A(:,i-listInd{j}(1)+1); 
                    end
                end
            end
            
            A = [A; Atemp];
            b = [b; cZ.b];
        end
    end

    % compute zonotope enclosure of the polynomial zonotopes
    ind = find(prod(ones(size(temp.expMat))-mod(temp.expMat,2),1) == 1);

    c = temp.c + 0.5 * sum(temp.G(:,ind),2);
    G = temp.G;
    G(:,ind) = 0.5*G(:,ind);
    
    % remove all-zero constraints
    ind = find(sum(abs(A),2) == 0);
    A(ind,:) = [];
    b(ind) = [];

    % construct resulting constrained zonotope
    res = conZonotope(c,G,A,b);

end