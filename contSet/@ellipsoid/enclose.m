function E = enclose(E,varargin)
% enclose - encloses an ellipsoid and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in E, x2 \in E2, a \in [0,1] }
%    where E2 = M*E + Eplus
%
% Syntax:
%    E = enclose(E,E2)
%    E = enclose(E,M,Eplus)
%
% Inputs:
%    E - ellipsoid object
%    M - matrix for the linear transformation
%    Eplus - ellipsoid object added to the linear transformation
%    E2 - ellipsoid object (M*Eplus)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E = ellipsoid(eye(2));
%    M = [2,1;1,1];
%    Eplus = ellipsoid([1,0;0,3],[1;-1]);
%    E_enc = enclose(E,M,Eplus);
%    
%    figure; hold on;
%    plot(E);
%    plot(Eplus,[1,2],'r');
%    plot(E_enc,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   04-July-2022 (VG, CORAerror and argument check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
if length(varargin) == 1
    E2 = varargin{1};
    inputArgsCheck({{E,'att','ellipsoid','scalar'}; ...
                    {E2,'att','ellipsoid','scalar'}});
    equalDimCheck(E,E2);

elseif length(varargin) == 2
    M = varargin{1};
    Eplus = varargin{2};
    inputArgsCheck({{E,'att','ellipsoid','scalar'}; ...
                    {Eplus,'att','ellipsoid','scalar'}; ...
                    {M,'att','numeric',{'size',[dim(E),dim(Eplus)]}}});
    E2 = M*Eplus;
    
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% compute enclosure using convex hull
E = convHull(E,E2);

% ------------------------------ END OF CODE ------------------------------
