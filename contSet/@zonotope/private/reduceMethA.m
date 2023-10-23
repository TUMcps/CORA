function [Zred]=reduceMethA(Z,order)
% reduceMethA - apply exhaustive search
%
% Syntax:
%    [Zred,t]=reduceMethA(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       11-September-2008
% Last update:   06-March-2009
%                27-June-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize Z_red
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGenerators(Z,order);

if ~isempty(Gred)

    %reorder generators
    Gcells=reorderingFilter(Gred);

    %pick generator with the best volume
    Gtemp=volumeFilter(Gcells,Z);
    Gpicked=Gtemp{1};

    %Build transformation matrix P
    P=Gpicked;
    
    % map generators
    Gtrans = pinv(P)*Gred;

    % box generators
    Gbox=diag(sum(abs(Gtrans),2));

    % transform generators back
    Gred = P*Gbox;
end

%build reduced zonotope
Zred.c = center;
Zred.G = [Gunred,Gred];


% ------------------------------ END OF CODE ------------------------------
