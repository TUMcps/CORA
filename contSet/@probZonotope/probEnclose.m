function [pZ] = probEnclose(pZ,Ar,Rtrans)
% probEnclose - Encloses a probabilistic zonotope pZ and its map A*pZ up to
% the m-sigma bound
%
% Syntax:  
%    [pZ] = probEnclose(pZ,A,m)
%
% Inputs:
%    pZ - probabilistic zonotope object
%    Ar - interval matrix of the linear map multiplied by the time increment
%    r
%
% Outputs:
%    pZ - probabilistic zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      28-August-2007
% Last update:  29-February-2008
%               03-September-2009
%               08-September-2009
% Last revision: ---

%------------- BEGIN CODE --------------

%TO DO: change computation as it is in the dissertation!

%get dimension
dim=length(pZ.Z(:,1));

%retrieve uncertain zonotope deltaZ
Zaux=(expm(Ar)-eye(dim))*zonotope(pZ)+Rtrans;
deltaZ=enclose(zonotope(zeros(dim,1)),Zaux);
pZ=pZ+deltaZ;

%change probabilistic generators if det(A)<1
if trace(Ar)<0
    if pZ.gauss
        pZ.cov=expm(Ar)*pZ.cov*expm(Ar)';
    else
        pZ.g=expm(Ar)*pZ.g;
    end
end

%------------- END OF CODE --------------