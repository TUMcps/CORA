function []=writeZonoGen(Z, Zred, nameF)
% write a Zonotope Z and its reduced version to a file
% Author: Anna Kopetzki
% Written: June-2016

%get center and Generators of original zonotope
Zmatrix=get(Z,'Z');
cen=Zmatrix(:,1);
Zgen=Zmatrix(:,2:end);

cenS = mat2str(cen);
ZgenS = mat2str(Zgen);


%get center and Generators of reduced zonotope
Zredmatrix=get(Zred,'Z');
cenRed=Zredmatrix(:,1);
Zredgen=Zredmatrix(:,2:end);

cenRedS = mat2str(cenRed);
ZredgenS = mat2str(Zredgen);


% write zonotopes to a file
fileID = fopen(nameF, 'w');
fprintf(fileID, '# Original zonotope\n');
fprintf(fileID, '# center: %s\n', cenS);
fprintf(fileID, '# generators: %s\n', ZgenS);
fprintf(fileID, '# Reducded zonotope\n');
fprintf(fileID, '# center: %s\n', cenRedS);
fprintf(fileID, '# generators: %s\n', ZredgenS);
