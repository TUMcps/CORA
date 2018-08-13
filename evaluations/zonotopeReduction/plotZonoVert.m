function []=plotZonoVert(Z, Zred, V, nameF)
% plot a Zonotope Z and its reduced version Zred and its generators
% Author: Anna Kopetzki (adapted from Mathias Althoff)
% Written: April-2016

%get center and Generators of original zonotope
Zmatrix=get(Z,'Z');
cen=Zmatrix(:,1);
Zgen=Zmatrix(:,2:end);
%get center and Generators of reduced zonotope
Zredmatrix=get(Zred,'Z');
cenRed=Zredmatrix(:,1);
Zredgen=Zredmatrix(:,2:end);

%plot zonotopes and Generators
%h = figure();
h=figure('Visible','off');
Vz = vertices(Z,[1,2]);
Vzred = vertices(Zred, [1,2]);
plot(Vz, 'b');
plot(Vzred, 'r');
plot(vertices(V), 'g');


quiver(cen(1)*ones(1,length(Zgen)),cen(2)*ones(1, length(Zgen)),Zgen(1,:),Zgen(2,:),0,'b');
quiver(cenRed(1)*ones(1,length(Zredgen)),cenRed(2)*ones(1, length(Zredgen)),Zredgen(1,:),Zredgen(2,:),0,'r');

c=clock;
nameFig=sprintf('%s_%f.png',nameF, c(6));
saveas(h,nameFig,'png');
