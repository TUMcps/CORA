function writeHessianTensorFile(J2dyn,J2con,path,fname,vars,infsupFlag,options)
% writeHessianTensorFile - generates an mFile that allows to compute the
%    hessian tensor
%
% Syntax:
%    writeHessianTensorFile(J2dyn,J2con,path,fname,vars,infsupFlag,options)
%
% Inputs:
%    J2dyn - Hessians of differential equation
%    J2con - Hessians of constraint equation
%    path - path for saving the file
%    fname - function name for the Hessian file
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
% See also: derivatives

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

[taylMod,opt] = aux_setOptions(options);

% squeeze dynamic and constraint part
Hdyn = aux_squeeze(J2dyn);
Hcon = aux_squeeze(J2con);

% open file
fid = fopen([path filesep fname '.m'],'w');

try

% function arguments depending on occurring variable types
[argsOut,argsIn,cellSymVars] = aux_args(vars);
fprintf(fid,'function %s = %s%s\n\n',argsOut,fname,argsIn);

% write call to optimized function to reduce the number of interval operations
[outDyn,indDyn] = aux_funOptimize(fid,Hdyn,vars.x,opt,'Dyn',argsIn);
[outCon,indCon] = aux_funOptimize(fid,Hcon,vars.y,opt,'Con',argsIn);

aux_writeTensors(fid,Hdyn,indDyn,vars.x,opt,taylMod,infsupFlag,'Hf');
aux_writeTensors(fid,Hcon,indCon,vars.y,opt,taylMod,infsupFlag,'Hg');

% add empty line and 'end'
fprintf(fid,'\nend');

% create optimized function to reduce the number of interval operations
if opt
    aux_writefunOptimize(fid,path,'funOptimizeDyn.m',outDyn,vars.x,cellSymVars);
    aux_writefunOptimize(fid,path,'funOptimizeCon.m',outCon,vars.y,cellSymVars);
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

% pre-processing
function [taylMod,opt] = aux_setOptions(options)
% init options for tensor generation

% default options
taylMod = false;
opt = false;

% read out options
if isfield(options,'lagrangeRem')
    taylMod = isfield(options.lagrangeRem,'method') ...
            && ~strcmp(options.lagrangeRem.method,'interval');
    opt = isfield(options.lagrangeRem,'simplify') ...
            && strcmp(options.lagrangeRem.simplify,'optimize');
end

end

function H = aux_squeeze(J2)

H = cell(size(J2,1),1);
for k=1:size(J2,1)
    H{k} = squeeze(J2(k,:,:));
end

end

function [argsOut,argsIn,cellSymVars] = aux_args(vars)
% generate output/input arguments and cell-array of used variables

% output arguments
if ~isempty(vars.y)
    argsOut = '[Hf,Hg]';
else
    argsOut = 'Hf';
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

function aux_writefunOptimize(fid,path,fname,outFun,vars,cellSymVars)

if ~isempty(vars)
    pathOpt = fullfile(path,fname);
    matlabFunction(outFun,'File',pathOpt,'Vars',cellSymVars);

    % fix missing 'end' etc.
    aux_fix_optimizedFuncs(pathOpt,fid);

    % delete funOptimizeDyn|Con file
    delete([path filesep fname]);
end

end

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

function [outFun,ind] = aux_funOptimize(fid,H,vars,opt,suffix,argsIn)

% ensure object consistency
outFun = [];
ind = repmat({struct('row',[],'col',[],'index',[])},size(H));

if ~opt || isempty(vars)
    return
end

% dynamic part
fprintf(fid,'out%s = funOptimize%s%s;\n',suffix,suffix,argsIn);

% store indices of nonempty entries
counter = 1;
for i = 1:length(H)
    [r,c] = find(H{i});
    if ~isempty(r)
        ind{i}.row = r;
        ind{i}.col = c;
        ind{i}.index = counter:counter + length(r)-1;
        counter = counter + length(r);
        for j = 1:length(r)
            outFun = [outFun;H{i}(r(j),c(j))]; 
        end
    end
end

end

function aux_writeTensors(fid,H,ind,vars,opt,taylMod,infsupFlag,varName)

% dynamic part
if isempty(vars)
    fprintf(fid,'\n\n%s = {};\n\n',varName);
    return
end

for k=1:length(H)
    % get matrix size
    [rows,cols] = size(H{k});
    sparseStr = sprintf('sparse(%i,%i)',rows,cols);
    if infsupFlag 
        str = sprintf('%s{%i} = interval(%s,%s);',varName,k,sparseStr,sparseStr);
    else
        str = sprintf('%s{%i} = %s;',varName,k,sparseStr);
    end
    % write in file if Hessian is used as Lagrange remainder
    fprintf(fid, '\n%s\n', str);
    % write rest of matrix
    if ~opt
        writeSparseMatrix(fid,H{k},...
            sprintf('%s{%i}',varName,k),infsupFlag && taylMod);
    else
        writeSparseMatrixOptimized(fid,ind{k},...
            sprintf('%s{%i}',varName,k),taylMod);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
