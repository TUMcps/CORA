function createRemainderFile(obj,J2dyn,J2con,path,name,options)
% createRemainderFile - generates an mFile that allows to compute the
% lagrange remainder
%
% Syntax:  
%    createRemainderFile(obj)
%
% Inputs:
%    obj - nonlinear system object
%    J2dyn - second order jacobian (dynamic part)
%    J2con - second order jacobian (constraint part)
%    path - path for file location
%    name - name of the corresponding nonlinear system
%    options - struct containing the algorithm options
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
% Written:      29-October-2007 
% Last update:  03-February-2011
%               05-September-2012
%               12-October-2015
%               25-January-2016
%               05-August-2016
%               05-November-2017
%               12-November-2017
%               03-December-2017
%               25-January-2018 (NK)
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

%create symbolic variables
[vars,vars_der] = symVariables(obj,'LRbrackets');

if isempty(J2con) % no constraint equations
    vars.z = [vars.x;vars.u]; % variables
    vars_der.z = [vars_der.x;vars_der.u]; % derivatives
else % with constraint equations
    vars.z = [vars.x;vars.y;vars.u]; % variables
    vars_der.z = [vars_der.x;vars_der.y;vars_der.u]; % derivatives 
end

% read out options
replace = 0;

if isfield(options,'replacements')
    replace = 1;
    if ~isempty(vars.p)
        rep = options.replacements(vars.x,vars.u,vars.p);
    else
        rep = options.replacements(vars.x,vars.u);
    end
end

% no parameters
if isempty(vars.p)
    for k=1:length(vars.x)
        lf(k,1)=0.5*vars_der.z.'*squeeze(J2dyn(k,:,:))*vars_der.z;
    end
    % Lagrange remainder of constraints
    for k=1:length(vars.y)
        lg(k,1)=0.5*vars_der.z.'*squeeze(J2con(k,:,:))*vars_der.z;
    end
    % set arguments
    if isempty(J2con) % no constraint equations
        if isempty(vars.T)
            arguments = '(x,u,dx,du)';
        else
            arguments = '(x,u,dx,du,T)';
        end
        returns = 'lf';
    else
        if isempty(vars.T)
            arguments = '(x,y,u,dx,dy,du)';
        else
            arguments = '(x,y,u,dx,dy,du,T)';
        end
        returns = '[lf, lg]';
    end
else
%parameters included
    for k=1:length(vars.x)
        lf(k,1)=0.5*vars_der.z.'*squeeze(J2dyn{1}(k,:,:))*vars_der.z;
        for iParam=1:length(vars.p)
            lf(k,1)=lf(k,1) + vars.p(iParam)*(0.5*vars_der.z.'*squeeze(J2dyn{iParam+1}(k,:,:))*vars_der.z);
        end
    end
    % set arguments
    if isempty(vars.T)
        arguments = '(x,u,dx,du,p)';
    else
        arguments = '(x,u,dx,du,p,T)';
    end
    returns = 'lf';
end


%create file
fid = fopen([path '/lagrangeRemainder_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function ',returns,'=lagrangeRemainder_',name,arguments]);
fprintf(fid, '%s\n\n', 'lf=interval();');
if ~isempty(J2con) % if constraints exist
    fprintf(fid, '%s\n\n', 'lg=interval();');
end


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
    
    % substite replacements in the dynamic equations
    lf = subs(lf,rep,r);
end


% dynamic part
for k=1:length(vars.x)
    if taylMod
        str=['lf(',num2str(k),',1)=interval(',char(lf(k,1)),');'];
    else
        str=['lf(',num2str(k),',1)=',char(lf(k,1)),';'];
    end
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
        
    %write in file
    fprintf(fid, '%s\n', str);
end

% constraint part
for k=1:length(vars.y)
    str=['lg(',num2str(k),',1)=',char(lg(k,1)),';'];
    %generate left and right brackets
    str=strrep(str,'L','(');
    str=strrep(str,'R',')');
        
    %write in file
    fprintf(fid, '%s\n', str);
end

%close file
fclose(fid);

%------------- END OF CODE --------------