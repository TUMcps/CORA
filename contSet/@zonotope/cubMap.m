function res = cubMap(Z,varargin)
% cubMap - computes an enclosure of the set corresponding to the cubic 
%    multiplication of a zonotope with a third-order tensor
%
% Description:
%    Calculates the following set:
%    { z = (x' T x) * x | x \in Z }
%
%    If three polyZonotopes are provided, the function calculates the set:
%    { z = (x1' T x2) * x3 | x1 \in Z1, x2 \in Z2, x3 \in Z3 }
%
% Syntax:
%    res = cubMap(Z,T)
%    res = cubMap(Z,T,ind)
%    res = cubMap(Z1,Z2,Z3,T)
%    res = cubMap(Z1,Z2,Z3,T,ind)
%
% Inputs:
%    Z,Z1,Z2,Z3 - zonotope objects
%    T - third-order tensor
%    ind - cell-array containing the non-zero indices of the tensor
%
% Outputs:
%    res - zonotope object representing the set of the cubic mapping
%
% Example: 
%    % cubic multiplication
%    Z = zonotope([1;-1],[1 3 -2 -1; 0 2 -1 1]);
%    
%    T{1,1} = rand(2); T{1,2} = rand(2);
%    T{2,1} = rand(2); T{2,2} = rand(2);
% 
%    Zcub = cubMap(Z,T);
%
%    figure;
%    subplot(1,2,1);
%    plot(Z,[1,2],'FaceColor','r');
%    subplot(1,2,2);
%    plot(Zcub,[1,2],'FaceColor','b');
%
%    % mixed cubic multiplication
%    Z2 = zonotope([1;-1],[-1 3 -2 0 3; -2 1 0 2 -1]);
%    Z3 = zonotope([1;-1],[-3 2 0; 2 1 -1]);
%
%    ZcubMixed = cubMap(Z,Z2,Z3,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadMap

% Authors:       Niklas Kochdumper
% Written:       17-August-2018
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
        
        % assign input arguments
        Z2 = varargin{1};
        Z3 = varargin{2};
        T = varargin{3};
        
        % parse optional input arguments
        if nargin == 5
            ind = varargin{4}; 
        else
            temp = 1:size(T,2);
            ind = repmat({temp},[size(T,1),1]);
        end

        % check input arguments
        inputArgsCheck({{Z,'att','zonotope'};
                        {Z2,'att','zonotope'};
                        {Z3,'att','zonotope'};
                        {T,'att','cell'};
                        {ind,'att','cell'}});
        
        % mixed cubic multiplication
        res = aux_cubMapMixed(Z,Z2,Z3,T,ind);
        
    elseif nargin == 2 || nargin == 3
        % res = cubMap(Z,T)
        % res = cubMap(Z,T,ind)
        
        % assign input arguments
        T = varargin{1};
        
        % parse optional input arguments
        if nargin > 2
            ind = varargin{2}; 
        else
            temp = 1:size(T,2);
            ind = repmat({temp},[size(T,1),1]);
        end 
        
        % check input arguments
        inputArgsCheck({{Z,'att','zonotope'};
                        {T,'att','cell'};
                        {ind,'att','cell'}});

        % cubic multiplication
        res = aux_cubMapSingle(Z,T,ind);       
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_cubMapSingle(Z,T,ind)
% calulates the following set:      { z = (x' T x) * x | x \in Z }

    % initialize variables
    n = length(ind);
    N = size(Z.G,2)+1;
    Zcub = zeros(n,N^3);

    % loop over all system dimensions
    for i = 1:length(ind)

       listQuad = repmat({zeros(N,N)},[1,N]);

       % loop over all quadratic matrices: \sum_k (x' T_k x) * x_k 
       for k = 1:length(ind{i})

           % quadratic evaluation
           quadMat = [Z.c,Z.G]' * T{i,ind{i}(k)} * [Z.c,Z.G];

           % add up all entries that correspond to identical factors
           temp = tril(quadMat,-1);
           quadMat = quadMat - temp;
           quadMat = quadMat + temp';

           % multiply with the zonotope generators of the corresponding
           % dimension
           for j = 1:N
               Zj = [Z.c,Z.G];
               listQuad{j} = listQuad{j} + quadMat * Zj(ind{i}(k),j);
           end
       end 

       % add up all entries that belong to identical factors
       for k = 2:N

           % loop over all quadratic matrix rows whos factors already appear in
           % one of the previous quadratic matrices
           for j = 1:k-1

               % loop over all row entries
               for h = j:N
                  if h <= k
                      listQuad{j}(h,k) = listQuad{j}(h,k) + listQuad{k}(j,h);
                  else
                      listQuad{j}(k,h) = listQuad{j}(k,h) + listQuad{k}(j,h);
                  end
               end
           end
       end

       % half the entries for purely quadratic factors
       temp = diag(listQuad{1});
       listQuad{1}(1,1) = listQuad{1}(1,1) + 0.5*sum(temp(2:end));

       for k = 2:N
          listQuad{1}(k,k) = 0.5 * listQuad{1}(k,k); 
       end

       % summerize all identical factors in one matrix
       counter = 1;

       for k = 1:N

           % loop over all matrix rows that contain unique factors
           for j = k:N
               m = N-j+1;         % number of elements in the row
               Zcub(i,counter : counter + m - 1) = listQuad{k}(j,j:end);
               counter = counter + m;
           end
       end
    end

    % concatenate the generator matrices
    Zcub = Zcub(:,1:counter-1);

    % construct the resulting zonotope
    res = zonotope(Zcub);
end

function res = aux_cubMapMixed(Z1,Z2,Z3,T,ind)
% calculates the following set:
% { z = (x1' T x2) * x3 | x1 \in pZ1, x2 \in pZ2, x3 \in pZ3 }

    % initialize variables
    n = length(ind);
    N1 = size(Z1.G,2)+1;
    N2 = size(Z2.G,2)+1;
    N3 = size(Z3.G,2)+1;
    Nq = N1*N2;

    Zcub = zeros(n,N1*N2*N3);

    % loop over all system dimensions
    for i = 1:length(ind)

       % loop over all quadratic matrices: \sum_k (x1' T_k x2) * x3_k 
       for k = 1:length(ind{i})

           % quadratic evaluation
           quadMat = [Z1.c,Z1.G]' * T{i,ind{i}(k)} * [Z2.c,Z2.G];
           quadVec = reshape(quadMat,1,[]);   

           % multiply with Z3
           for j = 1:N3
               Z3j = [Z3.c,Z3.G];
               Zcub(i,(j-1)*Nq + 1 : j*Nq) = Zcub(i,(j-1)*Nq + 1 : j*Nq) + ...
                                            quadVec * Z3j(ind{i}(k),j);
           end
       end 
    end

    % construct the resulting zonotope
    res = zonotope(Zcub);
end

% ------------------------------ END OF CODE ------------------------------
