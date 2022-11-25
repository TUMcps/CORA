function createHessianTensorFile_powSys(J2dyn,path,name,Ycomplex)
% createHessianTensorFile - generates an mFile that allows to compute the
% hessian tensor
%
% Syntax:  
%    createHessianTensorFile(obj,path)
%
% Inputs:
%    obj - nonlinear system object
%    path - path for saving the file
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
% Written:      17-May-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    
% %load second order jacobian
% J2dyn = obj.derivative.secondOrder.dyn;
% J2con = obj.derivative.secondOrder.con;


for k=1:length(J2dyn(:,1,1))
    Hdyn{k} = squeeze(J2dyn(k,:,:));
end

fid = fopen([path,'/hessianTensor_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function [Hf,Hg]=hessianTensor_',name,'(x,y,u)']);

%dynamic part
for k=1:length(Hdyn)
    %get matrix size
    [rows,cols] = size(Hdyn{k});
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    str=['Hf{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
    %write in file
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    writeSparseMatrix(Hdyn{k},['Hf{',num2str(k),'}'],fid);
    
    disp(['dynamic dim ',num2str(k)]);
end

%obtain constants
Y = abs(Ycomplex);
delta_load = angle(Ycomplex);

%number of constraints
load busIEEE14_genParam
nrOfBuses = length(Y);
nrOfGenerators = length(p.M);
nrOfLoadBuses = nrOfBuses - nrOfGenerators;


%generate symbolic constraint states
for i=1:2*nrOfBuses
    command=['y(',num2str(i),',1)=sym(''yL',num2str(i),'R'');'];
    eval(command);
end 

%create E, V, and Theta values
E = y(1:nrOfGenerators);
V(nrOfGenerators + 1 : nrOfGenerators + nrOfLoadBuses, 1) = y(nrOfGenerators + 1 : nrOfGenerators + nrOfLoadBuses, 1);
Theta(1:nrOfBuses, 1) = y(nrOfBuses + 1 : 2*nrOfBuses, 1);


%write replacements
%active power
for i = 1:nrOfBuses
    for n = 1:nrOfBuses
        C(i,n) = cos(delta_load(i,n) + Theta(n) - Theta(i));
    end
end

%reactive power
for i = 1:nrOfBuses
    for n = 1:nrOfBuses
        S(i,n) = sin(delta_load(i,n) + Theta(n) - Theta(i));
    end
end

% %maximum value of voltages
% %E
% for i = 1:nrOfGenerators
%     E_max(i,1) = max(E(i,1));
% end
% %V
% for i = 1:nrOfLoadBuses
%     V_max(i,1) = max(V(i,1));
% end

%write precomputation to file: active power 
for i = 1:nrOfBuses
    for n = 1:nrOfBuses
        str = ['C(',num2str(i),',',num2str(n),') = ',bracketSubs(char(C(i,n))),';'];
        %write in file
        fprintf(fid, '%s\n', str);
    end
end

%write precomputation to file: reactive power 
for i = 1:nrOfBuses
    for n = 1:nrOfBuses
        str = ['S(',num2str(i),',',num2str(n),') = ',bracketSubs(char(S(i,n))),';'];
        %write in file
        fprintf(fid, '%s\n', str);
    end
end

%write Emax and Vmax values
for i = 1:nrOfGenerators
    str = ['E_max(',num2str(i),',1) = inf(',bracketSubs(char(E(i,1))),');'];
    %write in file
    fprintf(fid, '%s\n', str);
end
for i = 1:nrOfLoadBuses
    str = ['V_max(',num2str(i),',1) = inf(',bracketSubs(char(V(i,1))),');'];
    %write in file
    fprintf(fid, '%s\n', str);
end


%constraint part
for k=1:length(Hcon)
    %get matrix size
    [rows,cols] = size(Hcon{k}); 
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    str=['Hg{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
    %write in file
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    writeAlgebraicHessian(k);
    
    disp(['constraint dim ',num2str(k)]);
end

%close file
fclose(fid);


function writeSparseMatrix(M,var,fid)


%write each row
[row,col] = find(M~=0);

for i=1:length(row)
    iRow = row(i);
    iCol = col(i);
    str=bracketSubs(char(M(iRow,iCol)));
    str=[var,'(',num2str(iRow),',',num2str(iCol),') = ',str,';'];
    %write in file
    fprintf(fid, '%s\n', str);
end

function writeAlgebraicHessian(k)

%sector x-x
for i = 1:nrOfGenerators
    for j = 1:nrOfGenerators
        str=[char(1/p.X_m(i)),'*E_max(',num2str(i),',1)*V_max(',num2str(i),',1)'
    end
end




function [str]=bracketSubs(str)

%generate left and right brackets
str=strrep(str,'L','(');
str=strrep(str,'R',')');

%------------- END OF CODE --------------