function [vars,vars_der] = symVariables(varargin)
% symVariables - generates symbolic variables of a continuous system 
%
% Syntax:  
%    [vars,vars_der] = symVariables(varargin)
%
% Inputs:
%    obj - contDynamics object
%    type - defines if 'LR' brackets should be used
%
% Outputs:
%    x - symbolic state variables
%    u - symbolic input variables
%    y - symbolic constraint variables
%    p - symbolic parameters
%    dx - symbolic state deviation form linearization point
%    du - symbolic input deviation form linearization point
%    dy - symbolic constraint deviation form linearization point
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
% Written:      18-January-2008
% Last update:  06-July-2017
%               05-November-2017
%               14-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==1
    obj=varargin{1};
    type=[];
elseif nargin==2
    obj=varargin{1};
    type=varargin{2};
end

if strcmp(type,'LRbrackets')
    %generate symbolic states
    if isprop(obj,'dim') && obj.dim>0
        for i=1:obj.dim
            command=['x(',num2str(i),',1)=sym(''xL',num2str(i),'R'');'];
            eval(command);
            command=['dx(',num2str(i),',1)=sym(''dxL',num2str(i),'R'');'];
            eval(command);
        end
    else
        x = [];
        dx = [];
    end

    %generate symbolic inputs
    if isprop(obj,'nrOfInputs') && obj.nrOfInputs>0
        for i=1:obj.nrOfInputs 
            command=['u(',num2str(i),',1)=sym(''uL',num2str(i),'R'');'];
            eval(command);
            command=['du(',num2str(i),',1)=sym(''duL',num2str(i),'R'');'];
            eval(command);  
        end  
    else
        u = [];
        du = [];
    end
    
    %generate symbolic constraint states
    if isprop(obj,'nrOfConstraints') && obj.nrOfConstraints>0
        for i=1:obj.nrOfConstraints
            command=['y(',num2str(i),',1)=sym(''yL',num2str(i),'R'');'];
            eval(command);
            command=['dy(',num2str(i),',1)=sym(''dyL',num2str(i),'R'');'];
            eval(command);
        end 
    else
        y = [];
        dy = [];
    end
    
    %generate symbolic parameters
    if isprop(obj,'nrOfParam') && obj.nrOfParam>0
        for i=1:obj.nrOfParam
            command=['p(',num2str(i),',1)=sym(''pL',num2str(i),'R'');'];
            eval(command);
        end 
    else
        p = [];
    end
    
else
    %generate symbolic states
    if isprop(obj,'dim') && obj.dim>0
        for i=1:obj.dim
            command=['x(',num2str(i),',1)=sym(''x',num2str(i),''');'];
            eval(command);
            command=['dx(',num2str(i),',1)=sym(''dx',num2str(i),''');'];
            eval(command);
        end
    else
        x = [];
        dx = [];
    end

    %generate symbolic inputs
    if isprop(obj,'nrOfInputs') && obj.nrOfInputs>0
        for i=1:obj.nrOfInputs
            command=['u(',num2str(i),',1)=sym(''u',num2str(i),''');'];
            eval(command);
            command=['du(',num2str(i),',1)=sym(''du',num2str(i),''');'];
            eval(command);
        end
    else
        u = [];
        du = [];
    end
    
    %generate symbolic constraint states
    if isprop(obj,'nrOfConstraints') && obj.nrOfConstraints>0
        for i=1:obj.nrOfConstraints
            command=['y(',num2str(i),',1)=sym(''y',num2str(i),''');'];
            eval(command);
            command=['dy(',num2str(i),',1)=sym(''dy',num2str(i),''');'];
            eval(command);
        end  
    else
        y = [];
        dy = [];
    end
    
    %generate symbolic parameters
    if isprop(obj,'nrOfParam') && obj.nrOfParam
        for i=1:obj.nrOfParam
            command=['p(',num2str(i),',1)=sym(''p',num2str(i),''');'];
            eval(command);
        end   
    else
        p = [];
    end  
end

% combine variables
vars.x = x;
vars.u = u;
vars.y = y;
vars.p = p;

vars_der.x = dx;
vars_der.u = du;
vars_der.y = dy;

%------------- END OF CODE --------------