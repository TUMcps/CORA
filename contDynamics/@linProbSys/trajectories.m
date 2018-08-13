function [obj]=trajectories(obj,Zinit)
% trajectories - plots sample trajectories of a linear interval system
% object, starting in the vertices of Zinit
%
% Syntax:  
%    [obj]=trajectories(obj,Zinit)
%
% Inputs:
%    obj - linear interval system object
%    Zinit - initial zonotope set
%
% Outputs:
%    obj - linear interval system object
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author:       Matthias Althoff
% Written:      23-January-2007 
% Last update:  15-June-2016
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from object structure
r=obj.r;
dim=obj.dim;
U=obj.U;
C=mid(U);
B=0.5*(supremum(U)-infimum(U));

%generate of system matrices (specification has to be edited by hand)
[A,A_min,A_max]=system_matrices(system,...
   [3,3,...
    3,3]);

%get vertices of initial zonotope
V=extreme(convert2polytope(Z_init));

%plot each trajectory
for i=1:length(A)
    %determine minimum eigenfrequency that determines the frequency of the
    %sine wave, the system will be stimulated with
    [Y,W]=eig(A{i});
    frequencies=abs(imag(diag(W)));
    min_freq=min(frequencies);
    write2file(A{i},B,C,min_freq);
    for iVertex=1:length(V(:,1))
        x0=V(iVertex,:)';
        clear dEqu.m
        [t,x] = ode45(@dEqu,[0,5],x0);
        hold on
        plot(x(:,4),x(:,5),'g-');
    end
end

%auxiliary function:
%write m-file of a sample system from the linear interval system
function write2file(A,B,C,min_freq)
    fid = fopen('dEqu.m','w');
    fprintf(fid, '%s\n\n', 'function dx=dEqu(t,x)');
    t_rand=pi/2*rand(1)
    %generate string
    str=['dx=',mat2str(A),'*x+',mat2str(C),'+'...
        ,mat2str(B),'*sin(',num2str(min_freq),'*(t+',num2str(t_rand),'));'];
    %write in file
    fprintf(fid, '%s\n', str);
    %close file
    fclose(fid);
end

%------------- END OF CODE --------------