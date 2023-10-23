function Z = intervalMultiplication(Z,I)
% intervalMultiplication - computes the multiplication of an interval with
%    a zonotope
%
% Syntax:
%    Z = intervalMultiplication(Z,I)
%
% Inputs:
%    Z - zonotope object 
%    I - interval object
%
% Outputs:
%    Z - zonotope after multiplication of an interval with a zonotope
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    I = interval([0 1; 1 0], [1 2; 2 1]);
%    IZ = intervalMultiplication(Z,I);
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(IZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       26-July-2016 
% Last update:   20-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%get center of interval matrix
T=center(I);
%get radius of interval matrix
S=rad(I);
%auxiliary value
Zabssum=sum(abs([Z.c,Z.G]),2);

%compute new zonotope
if ~any(T)
    % no empty generators if interval matrix is symmetric
    Z.c = 0*Z.c;
    Z.G = diag(S*Zabssum);
else
    Z.c = T*Z.c;
    Z.G = [T*Z.G,diag(S*Zabssum)];
end

% ------------------------------ END OF CODE ------------------------------
