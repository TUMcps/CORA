function create3rdOrderTensorFile(J3dyn,J3con,path,name,vars,options)
% create3rdOrderTensorFile - generates an mFile that allows to compute the
% 3rd order terms 
%
% Syntax:  
%    create3rdOrderTensorFile(obj,path)
%
% Inputs:
%    J3dyn - symbolic third-order tensor
%    J3con - symbolic thrid-order tensor (constraints)
%    path - path for saving the file
%    name - name of the nonlinear function to which the 3rd order tensor belongs
%    vars - structure containing the symbolic variables
%    options - structure containing the algorithm options
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      22-August-2012
% Last update:  08-March-2017
%               12-November-2017
%               03-December-2017
%               24-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

% read out options
taylMod = 0;
replace = 0;
parallel = 0;

if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')
    if ~ismember(options.lagrangeRem.method,{'taylorModel','zoo'})
       error('Wrong value for setting "options.lagrangeRem.method"!');
    end
    taylMod = 1;
end
if isfield(options,'replacements')
    replace = 1;
    if ~isempty(vars.p)
        rep = options.replacements(vars.x,vars.u,vars.p);
    else
        rep = options.replacements(vars.x,vars.u);
    end
end
if isfield(options,'tensorParallel') && options.tensorParallel == 1
    parallel = 1; 
end

%rearrange dynamic part
for k=1:length(J3dyn(:,1,1,1))
    for l=1:length(J3dyn(1,:,1,1))
        Tdyn{k,l} = squeeze(J3dyn(k,l,:,:));
    end
end

% rearrange constraint part
for k=1:length(J3con(:,1,1,1))
    for l=1:length(J3con(1,:,1,1))
        Tcon{k,l} = squeeze(J3con(k,l,:,:));
    end
end

% create the file
fid = fopen([path '/thirdOrderTensor_',name,'.m'],'w');

% function arguments depending on occurring variable types
if isempty(vars.y)      % no constraints
   strHead = ['function Tf = thirdOrderTensor_',name]; 
else                    % constraints
   strHead = ['function [Tf,Tg] = thirdOrderTensor_',name];
end

strHead = [strHead,'(x,u'];

if ~isempty(vars.p)     % parameter system
    strHead = [strHead,',p'];
end

if isempty(vars.T)      % time continious system
   strHead = [strHead,')']; 
else                    % time discrete system
   strHead = [strHead,',T)'];
end

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
if options.tensorOrder == 3
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
for k=1:length(Tdyn(:,1))
    
    if parallel
        caseStr = ['case ',num2str(k)];
        fprintf(fid, '\n %s \n', caseStr);
    end
    
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
                writeSparseMatrixTaylorModel(Tdyn{k,l},['C{i}{',num2str(l),'}'],fid);
            else
                str = sprintf(initStr,l);
                fprintf(fid, '\n\n %s\n\n', str);
                writeSparseMatrix(Tdyn{k,l},['C{i}{',num2str(l),'}'],fid);
            end
        else
            if taylMod
                str = sprintf(initStr,k,l);
                fprintf(fid, '\n\n %s\n\n', str);
                writeSparseMatrixTaylorModel(Tdyn{k,l},['Tf{',num2str(k),',',num2str(l),'}'],fid);
            else
                str = sprintf(initStr,k,l);
                fprintf(fid, '\n\n %s\n\n', str);
                writeSparseMatrix(Tdyn{k,l},['Tf{',num2str(k),',',num2str(l),'}'],fid);
            end
        end

        disp(['dynamic index ',num2str(k),',',num2str(l)]);
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
if ~isempty(vars.y)
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

            disp(['dynamic index ',num2str(k),',',num2str(l)]);
        end
    end
end

%close file
fclose(fid);


%------------- END OF CODE --------------