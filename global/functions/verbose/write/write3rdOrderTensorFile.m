function write3rdOrderTensorFile(J3dyn,J3con,path,fname,vars,infsupFlag,options)
% write3rdOrderTensorFile - generates an mFile that allows to compute the
%    third-order terms 
%
% Syntax:
%    write3rdOrderTensorFile(J3dyn,J3con,path,fname,vars,infsupFlag,options)
%
% Inputs:
%    J3dyn - symbolic third-order tensor
%    J3con - symbolic third-order tensor (constraints)
%    path - path for saving the file
%    fname - function name for the third-order tensor file
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
% See also: derivatives

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       22-August-2012
% Last update:   08-March-2017
%                12-November-2017
%                03-December-2017
%                24-January-2018 (NK)
%                13-March-2020 (NK, implemented options.simplify = optimize)
%                01-February-2021 (MW, add infsupFlag for different filenames)
% Last revision: 09-October-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% read out options
[taylMod,replace,parallel,opt,rep] = aux_setOptions(options,vars);

% rearrange dynamic and constraint part (4D -> 2D cell array of 2D)
Tdyn = aux_rearrange(J3dyn);
Tcon = aux_rearrange(J3con);
% sizes
[n_dyn,m_dyn] = size(Tdyn);

% open file
fid = fopen([path filesep fname '.m'],'w');

% try-catch block to ensure that file is closed
try

% function arguments depending on occurring variable types
[argsOut,argsIn,cellSymVars] = aux_args(vars);
fprintf(fid, 'function %s = %s%s\n\n',argsOut,fname,argsIn);

% write replacements
[rep,r] = aux_writeReplacements(fid,replace,rep);

% write call to optimized function to reduce the number of interval 
% arithmetic evaluations (e.g., sin(x(1)) is used multiple times)
[out,ind] = aux_funOptimize(Tdyn,opt);
if ~isempty(out)
    fprintf(fid, 'out = funOptimize%s;\n\n', argsIn);
end

% beginning of parallel execution
if parallel
    fprintf(fid, 'C = cell(%i,1);\n\n', n_dyn);
    fprintf(fid, 'parfor i = 1:%i\n\n', n_dyn);
    fprintf(fid, '%sswitch i\n', aux_tabs(1));
end

% pre-compute initialization string
initStr = aux_initStr(Tdyn,parallel,infsupFlag);

% dynamic part
indZero = false(n_dyn,m_dyn);

for k=1:n_dyn
    if parallel
        fprintf(fid,'\n%scase %i\n',aux_tabs(2),k);
    end
    
    for l=1:m_dyn
        % substitute all replacements
        if replace
            Tdyn{k,l} = subs(Tdyn{k,l},rep,r);
        end
        
        % write matrix
        if parallel
            fprintf(fid, '\n%s\n', sprintf(initStr,l));
            indZero(k,l) = writeSparseMatrix(fid,Tdyn{k,l},...
                sprintf('%sC{i}{%i}',aux_tabs(3),l),taylMod);
        elseif opt
            fprintf(fid, '\n%s\n', sprintf(initStr,k,l));
            indZero(k,l) = writeSparseMatrixOptimized(fid,ind{k,l},...
                sprintf('Tf{%i,%i}',k,l),taylMod);
        else
            fprintf(fid, '\n%s\n', sprintf(initStr,k,l));
            indZero(k,l) = writeSparseMatrix(fid,Tdyn{k,l},...
                sprintf('Tf{%i,%i}',k,l),taylMod);
        end

        if options.verbose
            fprintf('     .. dynamic index %i,%i\n',k,l);
        end
    end
end

% end of parallel execution (rewrite to Tf)
if parallel
    fprintf(fid, '\n%send\n', aux_tabs(1));
    fprintf(fid, 'end\n\n');
    fprintf(fid, 'Tf = cell(%i,%i);\n',n_dyn,m_dyn);
    fprintf(fid, 'for i=1:%i\n', n_dyn);
    fprintf(fid, '%sfor j=1:%i\n', aux_tabs(1), m_dyn);
    fprintf(fid, '%sTf{i,j} = C{i}{j};\n', aux_tabs(2));
    fprintf(fid, '%send\n', aux_tabs(1));
    fprintf(fid, 'end\n\n');
end

% constraint part
for k=1:size(Tcon,1)
    for l=1:size(Tcon,2)
        [rows,cols] = size(Tcon{k,l});
        sparseStr = sprintf('sparse(%i,%i)', rows, cols);
        str = sprintf('Tg{%i,%i} = interval(%s,%s);\n',k,l,sparseStr,sparseStr);
        fprintf(fid, '%s\n\n', str);
        writeSparseMatrix(fid,Tcon{k,l},sprintf('Tg{%i,%i}',k,l));

        if options.verbose
            fprintf('     .. dynamic index %i, %i\n',k,l);
        end
    end
end

% invert values to represent non-zero indices
indNonZero = ~indZero;
fprintf(fid,'\nind = cell(%i,1);\n',size(indNonZero,1));
for i=1:size(indNonZero,1)
    fprintf(fid,'ind{%i} = %s;\n',i,mat2str(find(indNonZero(i,:)')));
end

% properly end function
fprintf(fid,'\nend\n');

% create optimized function to reduce the number of interval operations
aux_writefunOptimize(fid,out,path,cellSymVars);

catch ME
    % close file
    fclose(fid);
    rethrow(ME);
end

% close file
fclose(fid);

end


% Auxiliary functions -----------------------------------------------------

% pre-processing
function [taylMod,replace,parallel,opt,rep] = aux_setOptions(options,vars)
% init options for tensor generation

% init as false
taylMod = false;
replace = false;
parallel = false;
opt = false;
rep = [];

% read out from field 'lagrangeRem'
if isfield(options,'lagrangeRem')
    taylMod = isfield(options.lagrangeRem,'method') ...
        && ~strcmp(options.lagrangeRem.method,'interval');
    if isfield(options.lagrangeRem,'replacements')
        replace = true;
        if ~isempty(vars.p)
            rep = options.lagrangeRem.replacements(vars.x,vars.u,vars.p);
        else
            rep = options.lagrangeRem.replacements(vars.x,vars.u);
        end
    end
    parallel = isfield(options.lagrangeRem,'tensorParallel') ...
        && options.lagrangeRem.tensorParallel == 1;
    opt = isfield(options.lagrangeRem,'simplify') ...
        && strcmp(options.lagrangeRem.simplify,'optimize');
end

end

function T = aux_rearrange(J3)
% rearrange 4D array into 2D cell array of matrices

T = cell(size(J3,1),size(J3,2));
for k=1:size(J3,1)
    for l=1:size(J3,2)
        T{k,l} = squeeze(J3(k,l,:,:));
    end
end

end


% function signature
function [argsOut,argsIn,cellSymVars] = aux_args(vars)
% generate the function signature of the generated file

% output arguments
if ~isempty(vars.y)
    argsOut = '[Tf,Tg,ind]';
else
    argsOut = '[Tf,ind]';
end

% input arguments and cell-array of symbolic variables
if ~isempty(vars.y)
    nonemptyVars = [true;true;true;false];
elseif ~isempty(vars.p)
    nonemptyVars = [true;false;true;true];
else
    nonemptyVars = [true;false;true;false];
end
argsIn = {'x','y','u','p'};
argsIn = ['(' strjoin(argsIn(nonemptyVars),',') ')'];

cellSymVars = {vars.x,vars.y,vars.u,vars.p};
cellSymVars = cellSymVars(nonemptyVars);

end

function initStr = aux_initStr(Tdyn,parallel,infsupFlag)
% initialization string

[rows,cols] = size(Tdyn{1,1});
sparseStr = sprintf('sparse(%i,%i)',rows,cols);

if infsupFlag
    if parallel
        initStr = sprintf('%sC{i}{%%i} = interval(%s,%s);',...
            aux_tabs(3),sparseStr,sparseStr);
    else
        initStr = sprintf('Tf{%%i,%%i} = interval(%s,%s);',sparseStr,sparseStr);
    end
else
    if parallel
        initStr = sprintf('%sC{i}{%%i} = %s;',aux_tabs(3),sparseStr);
    else
        initStr = sprintf('Tf{%%i,%%i} = %s;',sparseStr);
    end
end

end

function str = aux_tabs(n)
% get tabs (two whitespaces) for pretty indenting of generated file

str = repmat('  ',1,n);

end

% replacements
function [rep,r] = aux_writeReplacements(fid,replace,rep)

r = [];
if replace
    % generate symbolic variables for the replacements (hack because
    %   sym('rL%dR',[1,1])
    % yields rL1R1 instead of rL1R... which is annoying)
    r = sym('rL%dR',[length(rep)+1,1]);
    r = r(1:end-1);
    
    if size(rep,1) == 1
        rep = transpose(rep); 
    end
    
    % write the replacements to the file
    fprintf(fid, '%% replacements\n');
    writeMatrix(fid,rep,'r','BracketSubs',true);
end

end


% functions for simplify = 'optimize'
function [out,ind] = aux_funOptimize(Tdyn,opt)

% ensure object consistency
out = [];
ind = repmat({struct('row',[],'col',[],'index',[])},size(Tdyn));

if ~opt
    return
end
    
% store indices of nonempty entries
counter = 1;
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

end

function aux_writefunOptimize(fid,out,path,cellSymVars)

% out is always empty if simplify is not 'optimize' (opt == true)
if isempty(out)
    return
end
    
% create file with optimized evaluation
pathTemp = fullfile(path,'funOptimize.m');
matlabFunction(out,'File',pathTemp,'Vars',cellSymVars);

% we paste the content of the generated file into the tensor file;
% in some MATLAB versions, matlabFunction does not add the keyword
% 'end' at the end, so we check whether it's there and append it if not
text = fileread(pathTemp);
fprintf(fid,'\n%s\n',text);
if ~contains(text,'end')
    fprintf(fid,'\nend\n\n');
end

end

% ------------------------------ END OF CODE ------------------------------
