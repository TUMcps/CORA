function createJacobianFile(Jdyn,Jcon,Jp,path,name,vars)
% createJacobianFile - generates an mFile that allows to compute the
%    jacobian at a certain state and input
%
% Syntax:
%    createJacobianFile(Jdyn,Jcon,Jp,path,name,vars)
%
% Inputs:
%    Jdyn - Jacobian of dynamic equation
%    Jcon - Jacobian of constraint equation
%    Jp - Jacobian w.r.t. parameters
%    path - path where the function should be created
%    name - function name for the file computing the Jacobian
%    vars - struct containing the symbolic variables
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

% Authors:       Matthias Althoff
% Written:       21-August-2012
% Last update:   05-August-2016
%                05-November-2017
%                03-December-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% open file
fid = fopen([path filesep name '.m'],'w');
try

% system has no uncertain parameters
if isempty(Jp)
    % write first line
    if isempty(vars.y) % no constraints
        fprintf(fid, '%s\n\n', ['function [A,B]=',name,'(x,u)']);
    else % with constraints
        fprintf(fid, '%s\n\n', ['function [A,B,C,D,E,F]=',name,'(x,y,u)']);
    end
  
    % DYNAMIC MATRICES
    % write "A=["
    fprintf(fid, '%s', 'A=[');
    % write rest of matrix
    if ~isempty(Jdyn.x)
        writeMatrix(Jdyn.x,fid);
    else
        fprintf(fid, '%s', '];');
    end

    % write "B=["
    fprintf(fid, '%s', 'B=[');
    % write rest of matrix
    if ~isempty(Jdyn.u)
        writeMatrix(Jdyn.u,fid);
    else
        fprintf(fid, '%s', '];');
    end
    
    if ~isempty(vars.y)
        % write "C=["
        fprintf(fid, '%s', 'C=[');
        % write rest of matrix
        if isfield(Jdyn,'y') && ~isempty(Jdyn.y)
            writeMatrix(Jdyn.y,fid);
        else
            fprintf(fid, '%s', '];');
        end

        % INPUT MATRICES
        % write "D=["
        fprintf(fid, '%s', 'D=[');
        % write rest of matrix
        if ~isempty(Jcon) && isfield(Jcon,'x') && ~isempty(Jcon.x)
            writeMatrix(Jcon.x,fid);
        else
            fprintf(fid, '%s', '];');
        end

        % write "E=["
        fprintf(fid, '%s', 'E=[');
        % write rest of matrix
        if ~isempty(Jcon) && isfield(Jcon,'u') && ~isempty(Jcon.u)
            writeMatrix(Jcon.u,fid);
        else
            fprintf(fid, '%s', '];');
        end

        % write "F=["
        fprintf(fid, '%s', 'F=[');
        % write rest of matrix
        if ~isempty(Jcon) && isfield(Jcon,'y') && ~isempty(Jcon.y)
            writeMatrix(Jcon.y,fid);
        else
            fprintf(fid, '%s', '];');
        end
    end

% system has uncertain parameters
else
    % write first line
    fprintf(fid, '%s\n\n', ['function [A,B]=',name,'(x,u,p)']);
    
    % SYSTEM MATRICES
    for iMatrix = 1:length(Jp.x)
        % write "A{i}=["
        fprintf(fid, '%s', 'A{', num2str(iMatrix),'}=[');
        % write rest of matrix
        writeMatrix(Jp.x{iMatrix},fid);
    end


    % INPUT MATRICES
    for iMatrix = 1:length(Jp.u)
        % write "B{i}=["
        fprintf(fid, '%s', 'B{', num2str(iMatrix),'}=[');
        % write rest of matrix
        writeMatrix(Jp.u{iMatrix},fid);
    end
end

catch ME
    % close file
    fclose(fid);

    rethrow(ME)
end

% close file
fclose(fid);

% ------------------------------ END OF CODE ------------------------------
