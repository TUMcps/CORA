function createRemainderFile_certainInput(obj)
% createRemainderFile - generates an mFile that allows to compute the
% lagrange remainder
%
% Syntax:  
%    createRemainderFile(obj)
%
% Inputs:
%    obj - nonlinear system object
%
% Outputs:
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      29-October-2007 
% Last update:  03-February-2011
% Last revision:---

%------------- BEGIN CODE --------------

    
%load second order jacobian
J2=obj.jacobians.secondOrder;

[~,~,~,~,dx,du] = symVariables(obj,'LRbrackets');
du = zeros(length(du),1);

for k=1:length(J2(:,1,1))
    f(k,1)=0.5*[dx;du].'*squeeze(J2(k,:,:))*[dx;du];
end

%simplify expression
f=simple(f);

%use variable precision arithmetic!
%f=vpa(f);


fid = fopen([coraroot '/contDynamics/stateSpaceModels/lagrangeRemainder.m'],'w');
fprintf(fid, '%s\n\n', 'function f=lagrangeRemainder(x,u,dx,du)');
for k=1:length(J2(:,1,1))
    str=['f(',num2str(k),',1)=',char(f(k,1)),';'];
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
    %replace zeros by interval-zeros
    str=strrep(str,'=0;','=interval(0,0);');
    
    %for tank example only!!!
    for i=1:length(J2(:,1,1))
        %replace ^3/2 by sqrt(^3)
        str=strrep(str,['x(',num2str(i),')^(3/2)'],['sqrt(x(',num2str(i),')^3)']);
    end    
    
    %write in file
    fprintf(fid, '%s\n', str);
end

%close file
status = fclose(fid);

%------------- END OF CODE --------------