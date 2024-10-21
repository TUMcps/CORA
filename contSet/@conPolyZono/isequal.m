function res = isequal(cPZ,S,varargin)
% isequal - checks if a constrained polynomial zonotope is equal to another
%    set or point
%
% Syntax:
%    res = isequal(cPZ,S)
%    res = isequal(cPZ,S,tol)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object or numeric
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[1 0 1;0 -1 1],[1 0 2;0 1 1],[1 -0.5], ...
%                       -1,[0 2;1 0],[0.4 0;0.1 1]);
%    cPZ2 = conPolyZono([0;0],[1 1 0;1 0 -1],[2 1 0;1 0 1],[0.5 -1], ...
%                       1,[2 0;0 1],[0 0.4;1 0.1]);
%    isequal(cPZ1,cPZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isequal, zonotope/isequal

% Authors:       Niklas Kochdumper
% Written:       27-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{cPZ,'att',{'conPolyZono','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% ensure that numeric is second input argument
[cPZ,S] = reorderNumeric(cPZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cPZ.precedence
    res = isequal(S,cPZ,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(cPZ,S,true)
    res = false;
    return
end

% conPolyZono-conPolyZono case
if isa(S,'conPolyZono')
    res = aux_isequal_conPolyZono(cPZ,S,tol);
    return
end

throw(CORAerror('CORA:noops',cPZ,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_conPolyZono(cPZ,S,tol)

% initialize result
res = false;

% remove redundancies
cPZ = compact_(cPZ,'all',eps);
S = compact_(S,'all',eps);

% compare number of constraints (quick check)
if size(cPZ.A,1) ~= size(S.A)
    return;
end

% compare number of generators (quick check)
if size(cPZ.G,2) ~= size(S.G,2) || ... 
    size(cPZ.A,2) ~= size(S.A,2) || ...
    size(cPZ.GI,2) ~= size(S.GI,2)
    return 
end

% compare identifier vectors
temp1 = sort(cPZ.id); temp2 = sort(unique([cPZ.id;S.id]));
E1 = cPZ.E; E2 = S.E; 
EC1 = cPZ.EC; EC2 = S.EC;

if length(temp1) ~= length(temp2) || ~all(temp1 == temp2)
    return;
elseif ~all(cPZ.id == S.id)
    [~,E1,E2] = mergeExpMatrix(cPZ.id,S.id,cPZ.E,S.E);
    if ~isempty(cPZ.A)
        [~,EC1,EC2] = mergeExpMatrix(cPZ.id,S.id,cPZ.EC,S.EC);
    end
end

% jointly compare dependent generators and exponent matrices
if ~compareMatrices([cPZ.G;E1],[S.G;E2])
    return
end

% compare constraints
if ~isempty(cPZ.A)
   
    % compare constraint exponent matrices
    if ~(all(all(withinTol(EC1,EC2,tol))))
        return
    end
    
    % compare constraints (consider that f(x)=0 is identical to -f(x)=0)
    A1 = [cPZ.b,cPZ.A]; A2 = [S.b,S.A];
    
    for i = 1:size(A1,1)
      if ~(all(all(withinTol(A1(i,:),A2(i,:),tol))))
         if ~(all(all(withinTol(-A1(i,:),A2(i,:),tol))))
              return;
          end
       end
    end
end

% compare center and independent generators
Z1 = zonotope(cPZ.c,cPZ.GI);
Z2 = zonotope(S.c,S.GI);
res = isequal(Z1,Z2);

end

% ------------------------------ END OF CODE ------------------------------
