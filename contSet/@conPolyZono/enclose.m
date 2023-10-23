function cPZ = enclose(cPZ,varargin)
% enclose - encloses a constrained polynomial zonotope and its affine
%    transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in cPZ, x2 \in cPZ2, a \in [0,1] }
%    where cPZ2 = M*cPZ + cPZplus
% 
% Syntax:
%    cPZ = enclose(cPZ,cPZ2)
%    cPZ = enclose(cPZ,M,Splus)
%
% Inputs:
%    cPZ - conPolyZono object
%    cPZ2 - conPolyZono object
%    M - matrix for the linear transformation
%    cPZplus - set added to the linear transformation
%
% Outputs:
%    cPZ - conPolyZono object enclosing cPZ and cPZ2
%
% Example: 
%    c = [0;0];
%    G = [1 0;0 1];
%    E = [1 0;0 1];
%    A = [1 -1];
%    b = 0;
%    EC = [2 0;0 1];
%    cPZ1 = conPolyZono(c,G,E,A,b,EC);
% 
%    M = [1 2;-1 0]; b = [2;3];
%    cPZ2 = M*cPZ1 + b;
%   
%    cPZ = enclose(cPZ1,cPZ2);
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor',[0.6 0.6 0.6],'Splits',10);
%    plot(cPZ1,[1,2],'FaceColor','r','Splits',8);
%    plot(cPZ2,[1,2],'FaceColor','b','Splits',8);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/enclose, zonotope/enclose

% Authors:       Niklas Kochdumper
% Written:       25-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    cPZ2 = varargin{1};
elseif nargin == 3
    if ~isa(varargin{2},'zonotope') && ~isa(varargin{2},'interval')
        throw(CORAerror('CORA:wrongValue','third',"'zonotope' or 'interval'"));
    end
    cPZ2 = (varargin{1}*cPZ) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% check if exponent matrices are identical
if ~all(size(cPZ.id) == size(cPZ2.id)) ||  ...
   ~all(cPZ.id == cPZ2.id) ||  ...
   ~all(size(cPZ.E) == size(cPZ2.E)) || ...
   ~all(all(cPZ.E == cPZ2.E)) || ...
   ~all(size(cPZ.EC) == size(cPZ2.EC)) || ...
   ~all(all(cPZ.EC == cPZ2.EC)) || ...
   ~all(size(cPZ.A) == size(cPZ2.A)) || ~all(all(cPZ.A == cPZ2.A))
    throw(CORAerror('CORA:specialError','Constraint polynomial zonotopes are not compatible!'));
end

% compute set with the enclose method for polynomial zonotopes
pZ1 = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.GI,cPZ2.E,cPZ2.id);
pZ = enclose(pZ1,pZ2);

% compute exponent matrix of constraint system
EC = [cPZ.EC; zeros(1,size(cPZ.EC,2))];

% construct resulting constrained polynomial zonotope
cPZ = conPolyZono(pZ.c,pZ.G,pZ.E,cPZ.A,cPZ.b,EC, pZ.GI,pZ.id);

% ------------------------------ END OF CODE ------------------------------
