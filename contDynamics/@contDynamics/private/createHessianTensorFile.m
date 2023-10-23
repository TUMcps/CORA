function createHessianTensorFile(J2dyn,J2con,path,name,vars,infsupFlag,options)
% createHessianTensorFile - generates an mFile that allows to compute the
%    hessian tensor
%
% Syntax:
%    createHessianTensorFile(J2dyn,J2con,path,name,vars,infsupFlag,options)
%
% Inputs:
%    J2dyn - Hessians of differential equation
%    J2con - Hessians of constraint equation
%    path - path for saving the file
%    name - function name for the Hessian file
%    vars - symbolic variables
%    infsupFlag - true if interval arithmetic, otherwise false
%    options - additional information for tensors
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

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       21-August-2012
% Last update:   08-March-2017
%                05-November-2017
%                03-December-2017
%                13-March-2020 (NK, implemented options.simplify = optimize)
%                01-February-2021 (MW, different filename due to infsupFlag)
%                01-June-2022 (MW, optimize for constraint part)
%                09-June-2022 (MA, considered case that only constraint part is non-empty)
%                10-October-2023 (TL, fix missing end in optimized function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default options
taylMod = false;
opt = false;

% read out options
if isfield(options,'lagrangeRem')
    temp = options.lagrangeRem;
    if isfield(temp,'method') && ~strcmp(temp.method,'interval')
        taylMod = true;
    end
    if isfield(temp,'simplify') && strcmp(temp.simplify,'optimize')
        opt = true;
    end
end

% init
Hdyn = cell(size(J2dyn,1),1); %length(vars.x)
Hcon = cell(size(J2con,1),1); %length(vars.y)

% squeeze dynamic part
for k=1:size(J2dyn,1) %length(vars.x)
    Hdyn{k} = squeeze(J2dyn(k,:,:));
end
% squeeze constraint part
for k=1:size(J2con,1) %length(vars.y)
    Hcon{k} = squeeze(J2con(k,:,:));
end

% open file
fid = fopen([path filesep name '.m'],'w');
try

% function arguments depending on occurring variable types
if isempty(vars.p)
    if isempty(vars.y) % no constraints
        argsin = '(x,u)';
        argsout = 'Hf';
    else % with constraints
        argsin = '(x,y,u)';
        argsout = '[Hf,Hg]';
    end
else
    argsin = '(x,u,p)';
    argsout = 'Hf';
end
fprintf(fid, '%s\n\n', ['function ' argsout '=' name argsin]);

% write call to optimized function to reduce the number of interval operations
if opt
    % symVars only used for call to matlabfunction (different from input
    % arguments to hessian file!)
    symVars = {'vars.x','vars.y','vars.u','vars.p'};
    symVars = strjoin(symVars([~isempty(vars.x) ~isempty(vars.y) ...
        ~isempty(vars.u) ~isempty(vars.p)]),',');
    argsinOpt = ['(' strrep(symVars,'vars.','') ')'];

    % dynamic part
    if ~isempty(vars.x)
        str = ['outDyn = funOptimizeDyn',argsinOpt,';'];
        fprintf(fid, '\n\n %s\n\n', str);
        
        % store indices of nonempty entries
        counter = 1;
        indDyn = cell(size(Hdyn));
        outDyn = [];
        
        for i = 1:length(Hdyn)
            [r,c] = find(Hdyn{i});
            if ~isempty(r)
                indDyn{i}.row = r;
                indDyn{i}.col = c;
                indDyn{i}.index = counter:counter + length(r)-1;
                counter = counter + length(r);
                for j = 1:length(r)
                    outDyn = [outDyn;Hdyn{i}(r(j),c(j))]; 
                end
            end
        end
    end

    % constraint part
    if ~isempty(vars.y)
        str = ['outCon = funOptimizeCon',argsinOpt,';'];
        fprintf(fid, '\n\n %s\n\n', str);
        
        % store indices of nonempty entries
        counter = 1;
        indCon = cell(size(Hcon));
        outCon = [];
        
        for i = 1:length(Hcon)
            [r,c] = find(Hcon{i});
            if ~isempty(r)
                indCon{i}.row = r;
                indCon{i}.col = c;
                indCon{i}.index = counter:counter + length(r)-1;
                counter = counter + length(r);
                for j = 1:length(r)
                    outCon = [outCon;Hcon{i}(r(j),c(j))]; 
                end
            end
        end
    end
end

%dynamic part
if ~isempty(vars.x)
    for k=1:length(Hdyn)
        % get matrix size
        [rows,cols] = size(Hdyn{k});
        sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
        if infsupFlag 
            str=['Hf{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
        else
            str=['Hf{',num2str(k),'} = ',sparseStr,';'];
        end
        % write in file if Hessian is used as Lagrange remainder
        fprintf(fid, '\n\n %s\n\n', str);
        % write rest of matrix
        if ~opt
            if infsupFlag && taylMod
                writeSparseMatrixTaylorModel(Hdyn{k},['Hf{',num2str(k),'}'],fid);
            else
                writeSparseMatrix(Hdyn{k},['Hf{',num2str(k),'}'],fid);
            end
        elseif ~isempty(indDyn{k})
            writeSparseMatrixOptimized(indDyn{k},['Hf{',num2str(k),'}'],fid,taylMod);
        end
        
        disp(['     .. dynamics dim ',num2str(k)]);
    end
else
    fprintf(fid,'\n\nHf = {};\n\n');
end

%constraint part
if ~isempty(vars.y)
    for k=1:length(Hcon)
        % get matrix size
        [rows,cols] = size(Hcon{k});
        sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
        if infsupFlag 
            str=['Hg{',num2str(k),'} = interval(',sparseStr,',',sparseStr,');'];
        else
            str=['Hg{',num2str(k),'} = ',sparseStr,';'];
        end
        % write in file if Hessian is used as Lagrange remainder
        fprintf(fid, '\n\n %s\n\n', str);
        % write rest of matrix
        if ~opt
            if infsupFlag && taylMod
                writeSparseMatrixTaylorModel(Hcon{k},['Hg{',num2str(k),'}'],fid);
            else
                writeSparseMatrix(Hcon{k},['Hg{',num2str(k),'}'],fid);
            end
        elseif ~isempty(indCon{k})
            writeSparseMatrixOptimized(indCon{k},['Hg{',num2str(k),'}'],fid,taylMod);
        end
        
        disp(['     .. constraint dim ',num2str(k)]);
    end
end

% add empty line and 'end'
fprintf(fid,'\nend');

% create optimized function to reduce the number of interval operations
if opt
    % dynamic part
    if ~isempty(vars.x)
        pathOpt = fullfile(path,'funOptimizeDyn.m');
        eval(['matlabFunction(outDyn,''File'',pathOpt,''Vars'',{',symVars,'});']);

        % fix missing 'end' etc.
        aux_fix_optimizedFuncs(pathOpt,fid);

        % delete funOptimizeDyn file
        delete([path filesep 'funOptimizeDyn.m']);
    end

    % constraint part
    if ~isempty(vars.y)
        pathOpt = fullfile(path,'funOptimizeCon.m');
        eval(['matlabFunction(outCon,''File'',pathOpt,''Vars'',{',symVars,'});']);

        % fix missing 'end' etc.
        aux_fix_optimizedFuncs(pathOpt,fid);
    
        % delete funOptimizeCon file
        delete([path filesep 'funOptimizeCon.m']);
    end
end

catch ME
    % close file
    fclose(fid);

    rethrow(ME)
end

% close file
fclose(fid);

end


% Auxiliary functions -----------------------------------------------------

function aux_fix_optimizedFuncs(pathOpt,fid)
    
    % print text from pathOpt to hessian tensor file
    text = fileread(pathOpt);

    % find last non-empty line
    idx = 0;
    while strcmp(text(end-idx),compose('\n')) || strcmp(text(end-idx),compose('\r'))
        idx = idx + 1;
    end
    
    % check if end is missing
    if ~strcmp(text(end-idx-2:end-idx),'end')
        % print with missing end
        fprintf(fid,'\n\n%s\nend\n',text);
    else
        % print without extra end
        fprintf(fid,'\n\n%s\n',text);
    end

end

% ------------------------------ END OF CODE ------------------------------
