function [qZsplit] = splitOneGen(qZ,genNr)
% splitOneGen - Splits one generator factor of a quadZonotope
%
% Syntax:  
%    [qZsplit] = splitOneGen(qZ,genNr)
%
% Inputs:
%    qZ - quadZonotope object
%    genNr - generator number for the splitting
%
% Outputs:
%    qZsplit - cell array of split quadZonotopes
%
% Example: 
%    ---
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%obtain required indices
indG = genNr;
indGsquare = genNr;

%obtain index combination of Gquad
depGens = length(qZ.G(1,:));
C = combinator(depGens,2,'c');
ind1 = find(C(:,1) == genNr);
ind1_corresponding = C(ind1,2);
ind2 = find(C(:,2) == genNr);
ind2_corresponding = C(ind2,1);

%indices of Gquad
try
    indGquad = [ind1;ind2];
catch
    disp('error');
end
indGquad_corresponding = [ind1_corresponding;ind2_corresponding];

%compute centers of splitted parallelpiped
c1 = qZ.c - 0.5*qZ.G(:,indG) + 0.25*qZ.Gsquare(:,indGsquare);
c2 = qZ.c + 0.5*qZ.G(:,indG) + 0.25*qZ.Gsquare(:,indGsquare);

%adjust generators
G1 = qZ.G;
G1(:,indG) = -0.5*qZ.G(:,indG) + 0.5*qZ.Gsquare(:,indGsquare);
G1(:,indGquad_corresponding) = qZ.G(:,indGquad_corresponding) - 0.5*qZ.Gquad(:,indGquad);
G2 = qZ.G;
G2(:,indG) = 0.5*qZ.G(:,indG) + 0.5*qZ.Gsquare(:,indGsquare);
G2(:,indGquad_corresponding) = qZ.G(:,indGquad_corresponding) + 0.5*qZ.Gquad(:,indGquad);

Gsquare1 = qZ.Gsquare;
Gsquare1(:,indGsquare) = 0.25*qZ.Gsquare(:,indGsquare);
Gsquare2 = Gsquare1;

Gquad1 = qZ.Gquad;
Gquad1(:,indGquad) = -0.5*qZ.Gquad(:,indGquad);
Gquad2 = qZ.Gquad;
Gquad2(:,indGquad) = 0.5*qZ.Gquad(:,indGquad);

%generate splitted parallelpipeds
qZsplit{1} = quadZonotope(c1,G1,Gsquare1,Gquad1,qZ.Grest);
qZsplit{2} = quadZonotope(c2,G2,Gsquare2,Gquad2,qZ.Grest);


%------------- END OF CODE --------------