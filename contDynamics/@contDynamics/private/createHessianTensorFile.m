function createHessianTensorFile(J2dyn,J2con,path,name,vars,infsupFlag,options)
% createHessianTensorFile - generates an mFile that allows to compute the
% hessian tensor
%
% Syntax:  
%    createHessianTensorFile(J2dyn,J2con,path,name,vars,infsupFlag,options)
%
% Inputs:
%    J2dyn - Hessians of differential equation
%    J2con - Hessians of constraint equation
%    path - path for saving the file
%    name - name of the nonlinear function to which the hessian belongs
%    vars - symbolic variables
%    infsupFlag - true if interval arithmetic, otherwise false
%    options - additional information for tensors
%
% Outputs:
%    - 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  08-March-2017
%               05-November-2017
%               03-December-2017
%               13-March-2020 (NK, implemented options.simplify = optimize)
%               01-February-2021 (MW, different filename due to infsupFlag)
% Last revision:---

%------------- BEGIN CODE --------------

% read out options
taylMod = 0;
opt = 0;

if isfield(options,'lagrangeRem')
    temp = options.lagrangeRem;
    if isfield(temp,'method') && ~strcmp(temp.method,'interval')
        taylMod = 1;
    end
    if isfield(temp,'simplify') && strcmp(temp.simplify,'optimize')
        opt = 1; 
    end
end

% init
Hdyn = cell(length(vars.x),1);
Hcon = cell(length(vars.y),1);

% squeeze dynamic part
for k=1:length(vars.x)
    Hdyn{k} = squeeze(J2dyn(k,:,:));
end
% squeeze constraint part
for k=1:length(vars.y)
    Hcon{k} = squeeze(J2con(k,:,:));
end

% different filename depending on whether interval arithmetic is used
if infsupFlag
    hessianname = 'hessianTensorInt_';
else
    hessianname = 'hessianTensor_';
end

fid = fopen([path filesep hessianname name '.m'],'w');
% function arguments depending on occurring variable types
if isempty(vars.p)
    if isempty(vars.y) % no constraints
        fprintf(fid, '%s\n\n', ['function Hf=' hessianname name '(x,u)']);
        args = '(x,u)';
        symVars = 'vars.x,vars.u';
    else % with constraints
        fprintf(fid, '%s\n\n', ['function [Hf,Hg]=' hessianname name '(x,y,u)']);
        args = '(x,y,u)';
        symVars = 'vars.x,vars.y,vars.u';
    end
else
    fprintf(fid, '%s\n\n', ['function Hf=' hessianname name '(x,u,p)']);
    args = '(x,u,p)';
    symVars = 'vars.x,vars.u,vars.p';
end

% write call to optimized function to reduce the number of interval operations
if opt  
    % write function
    str = ['out = funOptimize',args,';'];
    fprintf(fid, '\n\n %s\n\n', str);
    
    % store indices of nonempty enries
    counter = 1;
    ind = cell(size(Hdyn));
    out = [];
    
    for i = 1:length(Hdyn)
       [r,c] = find(Hdyn{i});
       if ~isempty(r)
           ind{i}.row = r;
           ind{i}.col = c;
           ind{i}.index = counter:counter + length(r)-1;
           counter = counter + length(r);
           for j = 1:length(r)
              out = [out;Hdyn{i}(r(j),c(j))]; 
           end
       end
    end  
end

%dynamic part
for k=1:length(Hdyn)
    %get matrix size
    [rows,cols] = size(Hdyn{k});
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    if infsupFlag 
        str=['Hf{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
    else
        str=['Hf{',num2str(k),'} = ',sparseStr,';'];
    end
    %write in file if Hessian is used as Lagrange remainder
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    if ~opt
        if infsupFlag && taylMod
            writeSparseMatrixTaylorModel(Hdyn{k},['Hf{',num2str(k),'}'],fid);
        else
            writeSparseMatrix(Hdyn{k},['Hf{',num2str(k),'}'],fid);
        end
    elseif ~isempty(ind{k})
        writeSparseMatrixOptimized(ind{k},['Hf{',num2str(k),'}'],fid,taylMod);
    end
    
    disp(['dynamics dim ',num2str(k)]);
end

%constraint part
for k=1:length(Hcon)
    %get matrix size
    [rows,cols] = size(Hcon{k});
    sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
    if infsupFlag 
        str=['Hg{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
    else
        str=['Hg{',num2str(k),'} = ',sparseStr,';'];
    end
    %write in file if Hessian is used as Lagrange remainder
    fprintf(fid, '\n\n %s\n\n', str);
    % write rest of matrix
    writeSparseMatrix(Hcon{k},['Hg{',num2str(k),'}'],fid);
    
    disp(['constraint dim ',num2str(k)]);
end

% create optimized function to reduce the number of interval operations
if opt
    % create file with optimized evaluation
    pathTemp = fullfile(path,'funOptimize.m');
    str = ['matlabFunction(out,''File'',pathTemp,''Vars'',{',symVars,'});'];
    eval(str);
    
    % read in text from the file
    text = fileread(pathTemp);

    % print text from file
    fprintf(fid, '%s',text);
end

%close file
fclose(fid);



%------------- END OF CODE --------------