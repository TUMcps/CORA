function res = isequal(cPZ1,cPZ2,varargin)
% isequal - checks if two constrained polynomial zonotopes are equal
%
% Syntax:
%    res = isequal(cPZ1,cPZ2)
%    res = isequal(cPZ1,cPZ2,tol)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object
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

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{cPZ1,'att','conPolyZono'};
                {cPZ2,'att','conPolyZono'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% initialize result
res = false;

% remove redundancies
cPZ1 = compact_(cPZ1,'all',eps);
cPZ2 = compact_(cPZ2,'all',eps);

% compare number of constraints (quick check)
if size(cPZ1.A,1) ~= size(cPZ2.A)
    return;
end

% compare number of generators (quick check)
if size(cPZ1.G,2) ~= size(cPZ2.G,2) || ... 
    size(cPZ1.A,2) ~= size(cPZ2.A,2) || ...
    size(cPZ1.GI,2) ~= size(cPZ2.GI,2)
    return 
end

% compare identifier vectors
temp1 = sort(cPZ1.id); temp2 = sort(unique([cPZ1.id;cPZ2.id]));
E1 = cPZ1.E; E2 = cPZ2.E; 
EC1 = cPZ1.EC; EC2 = cPZ2.EC;

if length(temp1) ~= length(temp2) || ~all(temp1 == temp2)
    return;
elseif ~all(cPZ1.id == cPZ2.id)
    [~,E1,E2] = mergeExpMatrix(cPZ1.id,cPZ2.id,cPZ1.E,cPZ2.E);
    if ~isempty(cPZ1.A)
        [~,EC1,EC2] = mergeExpMatrix(cPZ1.id,cPZ2.id,cPZ1.EC,cPZ2.EC);
    end
end

% jointly compare dependent generators and exponent matrices
if ~compareMatrices([cPZ1.G;E1],[cPZ2.G;E2])
    return
end

% compare constraints
if ~isempty(cPZ1.A)
   
    % compare constraint exponent matrices
    if ~(all(all(withinTol(EC1,EC2,tol))))
        return
    end
    
    % compare constraints (consider that f(x)=0 is identical to -f(x)=0)
    A1 = [cPZ1.b,cPZ1.A]; A2 = [cPZ2.b,cPZ2.A];
    
    for i = 1:size(A1,1)
      if ~(all(all(withinTol(A1(i,:),A2(i,:),tol))))
         if ~(all(all(withinTol(-A1(i,:),A2(i,:),tol))))
              return;
          end
       end
    end
end

% compare center and independent generators
Z1 = zonotope(cPZ1.c,cPZ1.GI);
Z2 = zonotope(cPZ2.c,cPZ2.GI);
res = isequal(Z1,Z2);

% ------------------------------ END OF CODE ------------------------------
