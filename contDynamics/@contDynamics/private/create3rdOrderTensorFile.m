function create3rdOrderTensorFile(J3dyn,J3con,path,name,vars,infsupFlag,options)
% create3rdOrderTensorFile - generates an mFile that allows to compute the
%    third-order terms 
%
% Syntax:
%    create3rdOrderTensorFile(J3dyn,J3con,path,name,vars,infsupFlag,options)
%
% Inputs:
%    J3dyn - symbolic third-order tensor
%    J3con - symbolic third-order tensor (constraints)
%    path - path for saving the file
%    name - function name for the third-order tensor file
%    vars - structure containing the symbolic variables
%    infsupFlag - true if interval arithmetic, otherwise false
%    options - structure containing the algorithm options
%
% Outputs:
%    -
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       22-August-2012
% Last update:   08-March-2017
%                12-November-2017
%                03-December-2017
%                24-January-2018 (NK)
%                13-March-2020 (NK, implemented options.simplify = optimize)
%                01-February-2021 (MW, add infsupFlag for different filenames)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out options
taylMod = false;
replace = false;
parallel = false;
opt = false;

if isfield(options,'lagrangeRem')
    temp = options.lagrangeRem;
    if isfield(temp,'method') && ~strcmp(temp.method,'interval')
        taylMod = true;
    end
    if isfield(temp,'replacements')
        replace = true;
        if ~isempty(vars.p)
            rep = temp.replacements(vars.x,vars.u,vars.p);
        else
            rep = temp.replacements(vars.x,vars.u);
        end
    end
    if isfield(temp,'tensorParallel') && temp.tensorParallel == 1
        parallel = true; 
    end
    if isfield(temp,'simplify') && strcmp(temp.simplify,'optimize')
        opt = true; 
    end
end

%rearrange dynamic part
for k=1:length(J3dyn(:,1,1,1))
    for l=1:length(J3dyn(1,:,1,1))
        Tdyn{k,l} = squeeze(J3dyn(k,l,:,:));
    end
end

% rearrange constraint part
if ~isempty(J3con)
    for k=1:length(J3con(:,1,1,1))
        for l=1:length(J3con(1,:,1,1))
            Tcon{k,l} = squeeze(J3con(k,l,:,:));
        end
    end
else
    Tcon = [];
end

% open file
fid = fopen([path filesep name '.m'],'w');
try

% function arguments depending on occurring variable types
if isempty(vars.y) || isempty(Tcon)
    % no constraints
    strHead = ['function [Tf,ind] = ' name]; 
else
    % constraints
    strHead = ['function [Tf,Tg,ind] = ' name];
end

symVars = {'vars.x','vars.y','vars.u','vars.p'};
symVars = strjoin(symVars([~isempty(vars.x) ~isempty(vars.y) ...
    ~isempty(vars.u) ~isempty(vars.p)]),',');
argsIn = ['(' strrep(symVars,'vars.','') ')'];
strHead = [strHead argsIn];

fprintf(fid, '%s\n\n',strHead);


% write replacements
if replace
    % generate symbolic variables for the replacements
    for i = 1:length(rep)
       command=['r(',num2str(i),',1)=sym(''rL',num2str(i),'R'');'];
       eval(command);    
    end

    if size(rep,1) == 1
       rep = transpose(rep); 
    end
    
    % write the replacements
    str = 'r = [';
    fprintf(fid,'\n\n %s', str);
    writeMatrix(rep,fid);
end
if ~exist('r','var')
    r = [];
end

% write call to optimized function to reduce the number of interval operations
if opt
    
    % store indices of nonempty entries
    counter = 1;
    ind = cell(size(Tdyn));
    out = [];
    
    for i = 1:size(Tdyn,1)
        for j = 1:size(Tdyn,2)
            [r,c] = find(Tdyn{i,j});
            if ~isempty(r)
                ind{i,j}.row = r;
                ind{i,j}.col = c;
                ind{i,j}.index = counter:counter + length(r)-1;
                counter = counter + length(r);
                for k = 1:length(r)
                    out = [out;Tdyn{i,j}(r(k),c(k))]; 
                end
            end
       end
    end  
    
    % write function
    if ~isempty(out)
        index = find(strHead == '(');
        args = strHead(index(1):end);
        str = ['out = funOptimize',args,';'];
        fprintf(fid, '\n\n %s\n\n', str);
    end
end


% beginning of parallel execution
if parallel
    dim = length(Tdyn(:,1));
    str = ['C = cell(',num2str(dim),',1);'];
    fprintf(fid, '\n\n %s\n\n', str);
    str = ['parfor i = 1:',num2str(dim)];
    fprintf(fid, '\n\n %s\n\n', str);
    str = 'switch i';
    fprintf(fid, '\n\n %s\n\n', str);
end


% precompute initialization string
[rows,cols] = size(Tdyn{1,1});
sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
if infsupFlag
    if parallel
        initStr = ['C{i}{%i} = interval(',sparseStr,',',sparseStr,');'];
    else
        initStr = ['Tf{%i,%i} = interval(',sparseStr,',',sparseStr,');'];
    end
else
    if parallel
        initStr = ['C{i}{%i} = ',sparseStr,';'];
    else
        initStr = ['Tf{%i,%i} = ',sparseStr,';'];
    end
end


% dynamic part
indNonZero = cell(length(Tdyn(:,1)),1);

for k=1:length(Tdyn(:,1))
    
    if parallel
        caseStr = ['case ',num2str(k)];
        fprintf(fid, '\n %s \n', caseStr);
    end
    
    % initialize vector of non-zero indices
    indNonZero{k} = [];
    
    for l=1:length(Tdyn(1,:))
        
        % substitude all replacements
        if replace
            Tdyn{k,l} = subs(Tdyn{k,l},rep,r);
        end
        
        % write matrix
        if parallel
            if taylMod
                str = sprintf(initStr,l);
                fprintf(fid, '\n\n %s\n\n', str);
                empty = writeSparseMatrixTaylorModel(Tdyn{k,l},['C{i}{',num2str(l),'}'],fid);
            else
                str = sprintf(initStr,l);
                fprintf(fid, '\n\n %s\n\n', str);
                empty = writeSparseMatrix(Tdyn{k,l},['C{i}{',num2str(l),'}'],fid);
            end
        else
            if opt
                str = sprintf(initStr,k,l);
                fprintf(fid, '\n\n %s\n\n', str);
                
                if ~isempty(ind{k,l})
                    writeSparseMatrixOptimized(ind{k,l},['Tf{',num2str(k),',',num2str(l),'}'],fid,taylMod);
                    empty = false;
                else
                    empty = true;
                end
            else
                if taylMod
                    str = sprintf(initStr,k,l);
                    fprintf(fid, '\n\n %s\n\n', str);
                    empty = writeSparseMatrixTaylorModel(Tdyn{k,l},['Tf{',num2str(k),',',num2str(l),'}'],fid);
                else
                    str = sprintf(initStr,k,l);
                    fprintf(fid, '\n\n %s\n\n', str);
                    empty = writeSparseMatrix(Tdyn{k,l},['Tf{',num2str(k),',',num2str(l),'}'],fid);        
                end
            end
        end
        
        if ~empty
           indNonZero{k} = [indNonZero{k}; l]; 
        end

        disp(['     .. dynamic index ',num2str(k),',',num2str(l)]);
    end
end

% end of parallel execution
if parallel
    fprintf(fid, '\n\n %s\n', 'end');
    fprintf(fid, '%s\n\n', 'end');
    str = ['for i=1:',num2str(dim)];
    fprintf(fid, '%s\n', str);
    str = ['for j=1:',num2str(length(Tdyn(1,:)))];
    fprintf(fid, '%s\n', str);
    str = ['Tf{i,j} = C{i}{j};'];
    fprintf(fid, '%s\n', str);
    fprintf(fid, '%s\n', 'end');
    fprintf(fid, '%s\n', 'end');
end

% constraint part
if ~isempty(Tcon)
    for k=1:length(Tcon(:,1))
        for l=1:length(Tcon(1,:))
            %get matrix size
            [rows,cols] = size(Tcon{k,l});
            sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
            str=['Tg{',num2str(k),',',num2str(l),'} = interval(',sparseStr,',',sparseStr,');'];
            %str=['Tg{',num2str(k),',',num2str(l),'} = infsup(',sparseStr,',',sparseStr,');']; %for INTLAB
            %write in file
            fprintf(fid, '\n\n %s\n\n', str);
            % write rest of matrix
            writeSparseMatrix(Tcon{k,l},['Tg{',num2str(k),',',num2str(l),'}'],fid);

            disp(['     .. dynamic index ',num2str(k),',',num2str(l)]);
        end
    end
end

% indices of non-zero elements
fprintf(fid,'\n ind = cell(%i,1);',length(indNonZero));

for i = 1:length(indNonZero)
    if ~isempty(indNonZero{i})
        fprintf(fid,'\n ind{%i} = [',i);
        for j = 1:length(indNonZero{i})-1
            fprintf(fid,'%i;',indNonZero{i}(j));
        end
        fprintf(fid,'%i];\n\n',indNonZero{i}(end));
    else
        fprintf(fid,'\n ind{%i} = [];\n\n',i);
    end
end

% properly end function
fprintf(fid, 'end\n\n');

% create optimized function to reduce the number of interval operations
if opt && ~isempty(out)
    % create file with optimized evaluation
    pathTemp = fullfile(path,'funOptimize.m');
    str = ['matlabFunction(out,''File'',pathTemp,''Vars'',{',symVars,'});'];
    eval(str);
    
    % read in text from the file
    text = fileread(pathTemp);

    % print text from file
    fprintf(fid, '%s\n',text);

    if ~contains(text, 'end')
        % returned function text doesn't contain the 'end' keyword 
        % in some Matlab versions
        fprintf(fid, 'end\n\n');
    end
end

catch ME
    % close file
    fclose(fid);

    rethrow(ME);
end

% close file
fclose(fid);


% ------------------------------ END OF CODE ------------------------------
