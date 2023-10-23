function res = cubMap(cZ,varargin)
% cubMap - computes an enclosure of the set corresponding to the cubic 
%    multiplication of a constrained zonotope with a third-order tensor
%
% Description:
%    Calculates the following set:
%    { z = (x' T x) * x | x \in cZ }
%
%    If three conZonotopes are provided, the function calculates the set:
%    { z = (x1' T x2) * x3 | x1 \in cZ1, x2 \in cZ2, x3 \in cZ3 }
%
% Syntax:
%    res = cubMap(cZ,T)
%    res = cubMap(cZ,T,ind)
%    res = cubMap(cZ1,cZ2,cZ3,T)
%    res = cubMap(cZ1,cZ2,cZ3,T,ind)
%
% Inputs:
%    cZ,cZ1,cZ2,cZ3 - conZonotope objects
%    T - third-order tensor
%    ind - cell-array containing the non-zero indices of the tensor
%
% Outputs:
%    res - conZonotope object representing the set of the cubic mapping
%
% Example: 
%    % cubic multiplication
%    Z = [0 3 0 1 -1;0 0 2 1 1];
%    A = [1 0 1 -1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    
%    T{1,1} = [1 -1; 0 1]; T{1,2} = [1 0; 0 -1];
%    T{2,1} = [1 0; 0 -1]; T{2,2} = [-1 1; 0 -1];
%
%    cZcub = cubMap(cZ,T);
%
%    figure;
%    subplot(1,2,1);
%    plot(cZ,[1,2],'FaceColor','r');
%    subplot(1,2,2);
%    plot(cZcub,[1,2],'FaceColor','b','Template',50);
%
%    % mixed cubic multiplication
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ2 = conZonotope(Z,A,b);
%    Z = [1 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ3 = conZonotope(Z,A,b);
%
%    ZcubMixed = cubMap(cZ,cZ2,cZ3,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadMap, zonotope/cubMap

% Authors:       Niklas Kochdumper
% Written:       02-November-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check number of input arguments
    if nargin < 2
        throw(CORAerror('CORA:notEnoughInputArgs',2));
    elseif nargin > 5
        throw(CORAerror('CORA:tooManyInputArgs',5));
    end

    % cubic multiplication or mixed cubic multiplication
    if nargin == 4 || nargin == 5
        % res = cubMap(cZ1,cZ2,cZ3,T)
        % res = cubMap(cZ1,cZ2,cZ3,T,ind)
        
        % assign input arguments
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

        % check input arguments
        inputArgsCheck({{cZ,'att','conZonotope'};
                        {cZ2,'att','conZonotope'};
                        {cZ3,'att','conZonotope'};
                        {T,'att','cell'};
                        {ind,'att','cell'}});
        
        % mixed cubic multiplication
        res = aux_cubMapMixed(cZ,cZ2,cZ3,T,ind);
        
    elseif nargin == 2 || nargin == 3
        % res = cubMap(cZ,T)
        % res = cubMap(cZ,T,ind)
        
        % assign input argument
        T = varargin{1};
        
        % parse optional input arguments
        if nargin > 2
            ind = varargin{2}; 
        else
            temp = 1:size(T,2);
            ind = repmat({temp},[size(T,1),1]);
        end 
        
        % check input arguments
        inputArgsCheck({{cZ,'att',{'conZonotope'},{''}};
                        {T,'att',{'cell'},{''}};
                        {ind,'att',{'cell'},{''}}});

        % cubic multiplication
        res = aux_cubMapSingle(cZ,T,ind);       
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_cubMapSingle(cZ,T,ind)
% calculates the following set:      { z = (x' T x) * x | x \in cZ }

    if isempty(cZ.A) 

        pZ = cubMap(zonotope(cZ.c,cZ.G),T,ind);
        res = conZonotope(pZ);
        
    else
        
        % rescale constrained zonotope to reduce over-approximation
        if ~isempty(cZ.A)
            cZ = rescale(cZ);
        end
        
        % compute cubic map for polynomial zonotope to keep track of
        % dependencies
        pZ = polyZonotope(zonotope(cZ.c,cZ.G));
        pZ = cubMap(pZ,T,ind);

        % compute constraint matrix from constraint matrix of original
        % constrained zonotope
        index = find(sum(pZ.E,1) == 1);
        A = zeros(size(cZ.A,1),size(pZ.G,2));
        b = cZ.b;

        for i = 1:size(pZ.E,1)
            ind_ = find(pZ.id == i);
            if ~isempty(ind_)
                ind = find(pZ.E(ind_,index) == 1);
                if ~isempty(ind)
                   A(:,index(ind)) = cZ.A(:,i); 
                end
            end
        end

        % compute zonotope enclosure of the polynomial zonotopes
        ind = find(prod(ones(size(pZ.E))-mod(pZ.E,2),1) == 1);

        c = pZ.c + 0.5 * sum(pZ.G(:,ind),2);
        G = pZ.G;
        G(:,ind) = 0.5*G(:,ind);
        
        % remove all-zero constraints
        ind = find(sum(abs(A),2) == 0);
        A(ind,:) = [];
        b(ind) = [];

        % construct resulting constrained zonotope
        res = conZonotope(c,G,A,b);
    end
    
end

function res = aux_cubMapMixed(cZ1,cZ2,cZ3,T,ind)
% calculates the following set:
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
    m1 = size(cZ1.G,2); ind1 = 1:m1;
    m2 = size(cZ2.G,2); ind2 = m1+1:m1+m2;
    m3 = size(cZ3.G,2); ind3 = m1+m2+1:m1+m2+m3;
    
    pZ1 = polyZonotope(cZ1.c,cZ1.G,[],eye(m1),ind1');
    pZ2 = polyZonotope(cZ2.c,cZ2.G,[],eye(m2),ind2');
    pZ3 = polyZonotope(cZ3.c,cZ3.G,[],eye(m3),ind3');
    
    pZ = cubMap(pZ1,pZ2,pZ3,T,ind);

    % compute constraint matrix from con. mat. of original conZonotope
    index = find(sum(pZ.E,1) == 1);
    A = []; b = [];
    
    list = {cZ1,cZ2,cZ3};
    listInd = {ind1,ind2,ind3};
    
    for j = 1:length(list)
    
        cZ = list{j};
        
        if ~isempty(cZ.A)
            
            Atemp = zeros(size(cZ.A,1),size(pZ.G,2));

            for i = listInd{j}
                ind_ = find(pZ.id == i);
                if ~isempty(ind_)
                    ind = find(pZ.E(ind_,index) == 1);
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
    ind = find(prod(ones(size(pZ.E))-mod(pZ.E,2),1) == 1);

    c = pZ.c + 0.5 * sum(pZ.G(:,ind),2);
    G = pZ.G;
    G(:,ind) = 0.5*G(:,ind);
    
    % remove all-zero constraints
    ind = find(sum(abs(A),2) == 0);
    A(ind,:) = [];
    b(ind) = [];

    % construct resulting constrained zonotope
    res = conZonotope(c,G,A,b);

end

% ------------------------------ END OF CODE ------------------------------
