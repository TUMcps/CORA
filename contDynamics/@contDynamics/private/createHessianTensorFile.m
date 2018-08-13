function createHessianTensorFile(J2dyn,J2con,path,name,vars,infsupFlag,options)
% createHessianTensorFile - generates an mFile that allows to compute the
% hessian tensor
%
% Syntax:  
%    createHessianTensorFile(obj,path)
%
% Inputs:
%    obj - nonlinear system object
%    path - path for saving the file
%    name - name of the nonlinear function to which the hessian belongs
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
% Written:      21-August-2012
% Last update:  08-March-2017
%               05-November-2017
%               03-December-2017
% Last revision:---

%------------- BEGIN CODE --------------

% ckeck if the taylor models or zoo-objects are used to evaluate the
% remainder
taylMod = 0;

if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')
    if ~ismember(options.lagrangeRem.method,{'taylorModel','zoo'})
       error('Wrong value for setting "options.lagrangeRem.method"!');
    end
    taylMod = 1;
end

% init
Hdyn = cell(length(vars.x));
Hcon = cell(length(vars.y));

% squeeze dynamic part
for k=1:length(vars.x)
    Hdyn{k} = squeeze(J2dyn(k,:,:));
end
% squeeze constraint part
for k=1:length(vars.y)
    Hcon{k} = squeeze(J2con(k,:,:));
end

fid = fopen([path '/hessianTensor_',name,'.m'],'w');
% function arguments depending on occurring variable types
if isempty(vars.p)
    if isempty(vars.y) % no constraints
        if isempty(vars.T)
            fprintf(fid, '%s\n\n', ['function Hf=hessianTensor_',name,'(x,u)']);
        else
            fprintf(fid, '%s\n\n', ['function Hf=hessianTensor_',name,'(x,u,T)']);
        end
    else % with constraints
        if isempty(vars.T)
            fprintf(fid, '%s\n\n', ['function [Hf,Hg]=hessianTensor_',name,'(x,y,u)']);
        else
            fprintf(fid, '%s\n\n', ['function [Hf,Hg]=hessianTensor_',name,'(x,y,u,T)']);
        end
    end
else
    if isempty(vars.T)
        fprintf(fid, '%s\n\n', ['function Hf=hessianTensor_',name,'(x,u,p)']);
    else
        fprintf(fid, '%s\n\n', ['function Hf=hessianTensor_',name,'(x,u,p,T)']);
    end
end

%dynamic part
for k=1:length(Hdyn)
    %get matrix size
    [rows,cols] = size(Hdyn{k});
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    if infsupFlag 
        str=['Hf{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
        %str=['Hf{',num2str(k),'} = infsup(',sparseStr,',',sparseStr,');']; %for INTLAB
    else
        str=['Hf{',num2str(k),'} = ',sparseStr,';'];
    end
    %write in file if Hessian is used as Lagrange remainder
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    if infsupFlag && taylMod
        writeSparseMatrixTaylorModel(Hdyn{k},['Hf{',num2str(k),'}'],fid);
    else
        writeSparseMatrix(Hdyn{k},['Hf{',num2str(k),'}'],fid);
    end
    
    disp(['dynamic dim ',num2str(k)]);
end

%constraint part
for k=1:length(Hcon)
    %get matrix size
    [rows,cols] = size(Hcon{k});
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    if infsupFlag 
        str=['Hg{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
        %str=['Hg{',num2str(k),'} = infsup(',sparseStr,',',sparseStr,');']; %for INTLAB
    else
        str=['Hg{',num2str(k),'} = ',sparseStr,';'];
    end
    %write in file if Hessian is used as Lagrange remainder
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    writeSparseMatrix(Hcon{k},['Hg{',num2str(k),'}'],fid);
    
    disp(['dynamic dim ',num2str(k)]);
end


%close file
fclose(fid);



%------------- END OF CODE --------------