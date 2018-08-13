function [x_0,x_c,A,u,t_0,a_s,B_s,C_s] = symVariables(varargin)
% symVariables - generates symbolic variables of the switching surface
%
% Syntax:  
%    [t_0,a_s,B_s,C_s] = symVariables(varargin)
%
% Inputs:
%    obj - switching surface object
%    type - defines if 'LR' brackets should be used
%
% Outputs:
%    t_0 - switching time constant
%    a_s - linear switching time parameter
%    B_s - quadratic switching time parameter
%    C_s - cubic switching time parameter
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==1
    obj=varargin{1};
    type=[];
elseif nargin==2
    obj=varargin{1};
    type=varargin{2};
end

%switching time constant
syms t_0

if strcmp(type,'LRbrackets')
   
    %x0
    for i=1:obj.dim
        command=['x_0(',num2str(i),',1)=sym(''x_0L',num2str(i),'R'');'];
        eval(command);
    end
    
    %xc
    for i=1:obj.dim
        command=['x_c(',num2str(i),',1)=sym(''x_cL',num2str(i),'R'');'];
        eval(command);
    end
    
    %u
    for i=1:obj.dim
        command=['u(',num2str(i),',1)=sym(''uL',num2str(i),'R'');'];
        eval(command);
    end
    
    %A
    for i=1:obj.dim
        for j=1:obj.dim
            command=['A(',num2str(i),',',num2str(j),')=sym(''AL',num2str(i),'CO',num2str(j),'R'');'];
            eval(command);
        end
    end
    
    %linear switching time parameter
    for i=1:obj.dim
        command=['a_s(',num2str(i),',1)=sym(''a_sL',num2str(i),'R'');'];
        eval(command);
    end

    %quadratic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            command=['B_s(',num2str(i),',',num2str(j),')=sym(''B_sL',num2str(i),'CO',num2str(j),'R'');'];
            eval(command);
            %consider symmetric entries
            command=['B_s(',num2str(j),',',num2str(i),')=sym(''B_sL',num2str(i),'CO',num2str(j),'R'');'];
            eval(command);
        end
    end    
    
    %cubic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            for k=j:obj.dim
                command=['C_s(',num2str(i),',',num2str(j),',',num2str(k),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(i),',',num2str(k),',',num2str(j),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(j),',',num2str(i),',',num2str(k),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(j),',',num2str(k),',',num2str(i),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(k),',',num2str(j),',',num2str(i),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(k),',',num2str(i),',',num2str(j),')=sym(''C_sL',num2str(i),'CO',num2str(j),'CO',num2str(k),'R'');'];
                eval(command);
            end
        end
    end
 
    
else
    
    %x0
    for i=1:obj.dim
        command=['x_0(',num2str(i),',1)=sym(''x_0(',num2str(i),')'');'];
        eval(command);
    end
    
    %xc
    for i=1:obj.dim
        command=['x_c(',num2str(i),',1)=sym(''x_c(',num2str(i),')'');'];
        eval(command);
    end
    
    %u
    for i=1:obj.dim
        command=['u(',num2str(i),',1)=sym(''u(',num2str(i),')'');'];
        eval(command);
    end
    
    %A
    for i=1:obj.dim
        for j=1:obj.dim
            command=['A(',num2str(i),',',num2str(j),')=sym(''A(',num2str(i),'CO',num2str(j),')'');'];
            eval(command);
        end
    end
    
    %linear switching time parameter
    for i=1:obj.dim
        command=['a_s(',num2str(i),',1)=sym(''a_s(',num2str(i),')'');'];
        eval(command);
    end

    %quadratic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            command=['B_s(',num2str(i),'CO',num2str(j),',1)=sym(''B_s(',num2str(i),num2str(j),')'');'];
            eval(command);
            %consider symmetric entries
            command=['B_s(',num2str(j),'CO',num2str(i),',1)=sym(''B_s(',num2str(i),num2str(j),')'');'];
            eval(command);
        end
    end    
    
    %cubic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            for k=j:obj.dim
                command=['C_s(',num2str(i),'CO',num2str(j),'CO',num2str(k),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(i),'CO',num2str(k),'CO',num2str(j),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(j),'CO',num2str(i),'CO',num2str(k),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(j),'CO',num2str(k),'CO',num2str(i),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(k),'CO',num2str(j),'CO',num2str(i),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
                %consider symmetric entries
                command=['C_s(',num2str(k),'CO',num2str(i),'CO',num2str(j),',1)=sym(''C_s(',num2str(i),num2str(j),num2str(k),')'');'];
                eval(command);
            end
        end
    end
    
      
end


%------------- END OF CODE --------------