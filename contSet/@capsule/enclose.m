function C = enclose(C,varargin)
% enclose - encloses a capsule and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in C, x2 \in C2, a \in [0,1] }
%    where C2 = M*C + Cplus
%
% Syntax:
%    E = enclose(C,C2)
%    E = enclose(C,M,Cplus)
%
% Inputs:
%    C - capsule object
%    M - matrix for the linear transformation
%    Cplus - capsule object added to the linear transformation
%    C2 - capsule object (M*Cplus)
%
% Outputs:
%    C - capsule object
%
% Example: 
%    C = capsule([1;2],[1;1],0.4);
%    M = [2,1; 1,1];
%    Cplus = capsule([0;0],[0.1;0],0.1);
%    C_enc = enclose(C,M,Cplus);
%    
%    figure; hold on;
%    plot(C);
%    plot(Cplus,[1,2],'r');
%    plot(C_enc,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: capsule/convHull, ellipsoid/enclose

% Authors:       Mark Wetzlinger
% Written:       17-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
if length(varargin) == 1
    C2 = varargin{1};
    inputArgsCheck({{C,'att','capsule','scalar'}; ...
                    {C2,'att','capsule','scalar'}});
    equalDimCheck(C,C2);

elseif length(varargin) == 2
    M = varargin{1};
    Cplus = varargin{2};
    inputArgsCheck({{C,'att','capsule','scalar'}; ...
                    {Cplus,'att','capsule','scalar'}; ...
                    {M,'att','numeric',{'size',[dim(C),dim(Cplus)]}}});
    C2 = M*Cplus;
    
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% compute enclosure using convex hull
C = convHull(C,C2);

% ------------------------------ END OF CODE ------------------------------
