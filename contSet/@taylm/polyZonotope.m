function res = polyZonotope(obj)
% polyZonotope - convert a Taylor model to a polynomial zonotope
%
% Syntax:
%    res = polyZonotope(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - polyZonotope object
%
% Examples: 
%    % create a taylor model 
%    syms x y
%    func = [sin(x+y) + exp(-x) + x*y;cos(x*y) + y - x];
%    t = taylm(func,interval([1;3],[2;4]),4);
%
%    % convert to a polynomial zonotope and plot the result
%    pZ = polyZonotope(t);
%    plot(pZ,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadZonotope/polyZonotope, zonotope/polyZonotope

% Authors:       Niklas Kochdumper
% Written:       12-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


    % get names of all variables
    names = {};
    
    for i = 1:size(obj,1)
       temp = obj(i,1);
       names = [names, temp.names_of_var];
    end
    
    names = unique(names);

    % check if the taylor model is valid
    if size(obj,2) > 1
        throw(CORAerror('CORA:specialError',...
            'Invalid taylor model object that can not be converted!'));
    end
    
    % initialize polynomial zonotope matrices
    n = size(obj,1);
    E = [];
    G = [];
    GI = zeros(n,n);
    c = zeros(n,1);
    
    % loop over all taylor model dimensions
    for i = 1:n
        
       temp = obj(i,1);
       
       % determine indices of variables
       ind = zeros(length(temp.names_of_var),1);
       
       for j = 1:length(temp.names_of_var)
           for k = 1:length(names)
              if strcmp(names{k},temp.names_of_var{j})
                  ind(j) = k;
                  break;
              end
           end
       end
       
       % convert polynomial part
       coeff = temp.coefficients';
       e = temp.monomials(:,2:end)';
       
       if i == 1
          G = [G, [coeff;zeros(n-1,length(coeff))]];
       elseif i == n
          G = [G, [zeros(n-1,length(coeff));coeff]];
       else
          G = [G, [zeros(i-1,length(coeff));coeff;zeros(n-i,length(coeff))]]; 
       end
       
       Etemp = zeros(length(names),size(e,2));
       Etemp(ind,:) = e;
       E = [E,Etemp];       
       
       % convert remainder part
       c(i) = center(temp.remainder);
       GI(i,i) = rad(temp.remainder);
    end
    
    % add all constant monomials to the center vector
    temp = sum(E,1);
    ind = find(temp == 0);
    c = c + sum(G(:,ind),2);
    E(:,ind) = [];
    G(:,ind) = [];
    
    % construct resulting polynomial zonotope object
    res = polyZonotope(c,G,GI,E) + polyZonotope(zeros(n,1),[],[]);
    
end

% ------------------------------ END OF CODE ------------------------------
