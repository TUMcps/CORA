function res = cubMap(pZ,varargin)
% cubMap - computes the set corresponding to the cubic multiplication of a 
%          polynomial zonotope with a third-order tensor
%
% Description:
%   Calulates the following set:
%   { z = (x' T x) * x | x \in pZ }
%
%   If three polyZonotopes are provided, the function calculates the set:
%   { z = (x1' T x2) * x3 | x1 \in pZ1, x2 \in pZ2, x3 \in pZ3 }
%
% Syntax:  
%    res = cubMap(pZ,T)
%    res = cubMap(pZ1,pZ2,pZ3,T)
%    res = cubMap(pZ,T,ind)     
%    res = cubMap(pZ1,pZ2,pZ3,T,ind)
%
% Inputs:
%    pZ,pZ1,pZ2,pZ3 - polyŹonotope objects
%    T - third-order tensor
%    ind - cell-array containing the non-zero indizes of the tensor
%
% Outputs:
%    res - polyZonotope object representing the set of the cubic mapping
%
% Example: 
%    % cubic multiplication
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    
%    T{1,1} = [1 2; -1 2];
%    T{1,2} = [-3 0; 1 1];
%    T{2,1} = [2 0; -2 1];
%    T{2,2} = [-3 0; -21 -1];
%
%    pZcub = cubMap(pZ,T);
%
%    figure 
%    subplot(1,2,1)
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    subplot(1,2,2)
%    plot(pZcub,[1,2],'b','Filled',true,'EdgeColor','none');
%
%    % mixed cubic multiplication
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZcubMixed = cubMap(pZ,pZ2,pZ2,T);
%
%    figure
%    subplot(1,3,1);
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    subplot(1,3,2);
%    plot(pZ2,[1,2],'b','Filled',true,'EdgeColor','none');
%    subplot(1,3,3);
%    plot(pZcubMixed,[1,2],'Filled',true,'g','EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadMap, zonotope/cubMap

% Author:       Niklas Kochdumper
% Written:      17-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % cubic multiplication or mixed cubic multiplication
    if isa(varargin{1},'polyZonotope')
        
        % check user input
        if nargin < 4
           error('Wrong syntax for function cubMap!'); 
        end
        
        pZ2 = varargin{1};
        pZ3 = varargin{2};
        T = varargin{3};
        
        % parse optional input arguments
        if nargin > 4
           ind = varargin{4}; 
        else
           temp = 1:size(T,2);
           ind = repmat({temp},[size(T,1),1]);
        end 
        
        % mixed cubic multiplication
        res = cubMapMixed(pZ,pZ2,pZ3,T,ind);
        
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
        res = cubMapSingle(pZ,T,ind);       
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = cubMapSingle(pZ,T,ind)
% calulates the following set:      { z = (x' T x) * x | x \in pZ }

    % split into a zonotope Z that overapproximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    pZtemp = pZ;
    pZtemp.Grest = [];

    m = size(pZ.Grest,2);
    temp = max(pZ.id);
    id_ = (temp:temp+m-1)';

    Zrem = polyZonotope(0*pZ.c, pZ.Grest, [], eye(m), id_);

    % construct extended generator and exponent matrix (extended by center)
    Gext = [pZ.c, pZ.G];
    Eext = [zeros(size(pZ.expMat,1),1), pZ.expMat];

    % initialize the resulting generator and exponent matrix 
    N = size(Gext,2);
    n = length(ind);
    M = N*(N+1)/2;

    Equad = zeros(size(pZ.expMat,1),M);
    Ecub = zeros(size(Equad,1),size(Equad,2)*N);
    Gcub = zeros(n,size(Equad,2)*N);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N
        Equad(:,counter:counter+N-j) = Eext(:,j:N) + Eext(:,j)*ones(1,N-j+1);
        counter = counter + N - j + 1;
    end
    
    % create the exponent matrix that corresponds to the cubic map
    for j = 1:N
        Ecub(:,(j-1)*M + 1 : j*M) = Equad + Eext(:,j)*ones(1,M);
    end

    % loop over all dimensions
    for i = 1:length(ind)
        
       % initialize quadratic matrix
       Gquad = zeros(1,M);
        
       % loop over all quadratic matrices: \sum_k (pZ' T_k pZ) * pZ_k 
       for k = 1:length(ind{i})
           
           % quadratic evaluation
           quadMat = Gext'*T{i,ind{i}(k)}*Gext;

           % copy the result into the generator matrix
           counter = 1;

           for j = 1:N
               Gquad(counter:counter+N-j) = ...
                          [quadMat(j,j) , quadMat(j,j+1:N) + quadMat(j+1:N,j)'];
                      
               counter = counter + N - j + 1;
           end
           
           % cubic generator matrix: loop over all generators
           for j = 1:N
               Gcub(i,(j-1)*M + 1 : j*M) = Gcub(i,(j-1)*M + 1 : j*M) + ...
                                           Gquad * Gext(ind{i}(k),j);
           end
       end
    end

    % add up all generators that belong to identical exponents
    [ExpNew,Gnew] = removeRedundantExponents(Ecub,Gcub);
    
    % mixed multiplication with the zonotope from the independent terms
    if ~isempty(pZ.Grest) && ~all(all(pZ.Grest==0))
        
        % cubic multiplication
        pZrest = cubMap(Zrem,T,ind);
        
        pZrestList = cell(6,1);
        pZrestList{1} = cubMap(Zrem,Zrem,pZtemp,T,ind);
        pZrestList{2} = cubMap(Zrem,pZtemp,pZtemp,T,ind);
        pZrestList{3} = cubMap(Zrem,pZtemp,Zrem,T,ind); 
        pZrestList{4} = cubMap(pZtemp,Zrem,pZtemp,T,ind);
        pZrestList{5} = cubMap(pZtemp,Zrem,Zrem,T,ind);
        pZrestList{6} = cubMap(pZtemp,pZtemp,Zrem,T,ind);
               
        % sum of all sets
        pZsum = sum(pZrest,pZrestList);
        
        % zonotope over-approximation
        Zsum = zonotope(pZsum);
        
        Grest = generators(Zsum);
        c = center(Zsum);
        
    else
        Grest = [];
        c = zeros(n,1);
    end
    
    % construct the resulting polynomial zonotope
    res = polyZonotope(Gnew(:,1) + c, Gnew(:,2:end), Grest, ...
                       ExpNew(:,2:end), pZ.id);
end

function res = cubMapMixed(pZ1,pZ2,pZ3,T,ind)
% calculates the following set:
% { z = (x1' T x2) * x3 | x1 \in pZ1, x2 \in pZ2, x3 \in pZ3 }

    % bring the exponent matrices to a common representation
    [id_,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);
    [id,expMat1,expMat3] = mergeExpMatrix(id_,pZ3.id,expMat1,pZ3.expMat);
    [~,expMat2,~] = mergeExpMatrix(id_,pZ3.id,expMat2,pZ3.expMat);
    
    % split into a zonotope Z that over-approximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    pZtemp = pZ1;
    pZtemp.Grest = [];
    Z1 = zonotope(pZtemp);
    Zrem1 = zonotope([0*pZ1.c, pZ1.Grest]);
    
    pZtemp = pZ2;
    pZtemp.Grest = [];
    Z2 = zonotope(pZtemp);
    Zrem2 = zonotope([0*pZ2.c, pZ2.Grest]);
    
    pZtemp = pZ3;
    pZtemp.Grest = [];
    Z3 = zonotope(pZtemp);
    Zrem3 = zonotope([0*pZ3.c, pZ3.Grest]);

    % construct extended generator and exponent matrix (extended by center)
    Gext1 = [pZ1.c, pZ1.G];
    Eext1 = [zeros(size(expMat1,1),1), expMat1];
    
    Gext2 = [pZ2.c, pZ2.G];
    Eext2 = [zeros(size(expMat2,1),1), expMat2];
    
    Gext3 = [pZ3.c, pZ3.G];
    Eext3 = [zeros(size(expMat3,1),1), expMat3];

    % initialize the resulting generator and exponent matrix 
    N1 = size(Gext1,2);
    N2 = size(Gext2,2);
    N3 = size(Gext3,2);
    
    n = length(ind);
    M = N1*N2;

    Equad = zeros(size(expMat1,1),M);
    Ecub = zeros(size(Equad,1),size(Equad,2)*N3);
    Gcub = zeros(n,size(Equad,2)*N3);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N2
        Equad(:,counter:counter+N1-1) = Eext1 + Eext2(:,j)*ones(1,N1);
        counter = counter + N1 ;
    end
    
    % create the exponent matrix that corresponds to the cubic map
    for j = 1:N3
        Ecub(:,(j-1)*M + 1 : j*M) = Equad + Eext3(:,j)*ones(1,M);
    end

    % loop over all dimensions
    for i = 1:length(ind)
        
       % loop over all quadratic matrices: \sum_k (pZ' T_k pZ) * pZ_k 
       for k = 1:length(ind{i})
           
           % quadratic evaluation
           quadMat = Gext1'*T{i,ind{i}(k)}*Gext2;
           quadVec = reshape(quadMat,1,[]);

           
           % cubic generator matrix: loop over all generators
           for j = 1:N3
               Gcub(i,(j-1)*M + 1 : j*M) = Gcub(i,(j-1)*M + 1 : j*M) + ...
                                           quadVec * Gext3(ind{i}(k),j);
           end
       end
    end

    % add up all generators that belong to identical exponents
    [ExpNew,Gnew] = removeRedundantExponents(Ecub,Gcub);
    
    % mixed multiplication with the zonotope from the independent terms       
    zonoTemp = cubMap(Z1,Zrem2,Z3,T,ind) + ...
               cubMap(Zrem1,Z2,Z3,T,ind) + ...
               cubMap(Zrem1,Zrem2,Z3,T,ind) + ...
               cubMap(Z1,Z2,Zrem3,T,ind) + ... 
               cubMap(Z1,Zrem2,Zrem3,T,ind) + ...
               cubMap(Zrem1,Z2,Zrem3,T,ind) + ...
               cubMap(Zrem1,Zrem2,Zrem3,T,ind);

    Grest = generators(zonoTemp);
    c = center(zonoTemp);
    
    % remove zero generators
    Grest(:,sum(abs(Grest),1) == 0) = [];
    
    % construct the resulting polynomial zonotope
    res = polyZonotope(Gnew(:,1) + c, Gnew(:,2:end), Grest, ExpNew(:,2:end),id);

end

%------------- END OF CODE --------------