function res = test_requiredToolboxes
% test_requiredToolboxes - checks if required toolboxes are installed
%
% Syntax:  
%    res = test_requiredToolboxes
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      15-September-2016
% Last update:  04-May-2018
% Last revision:---


%------------- BEGIN CODE --------------

    % check if symbolic math toolbox is available
    res_partial(1) = license('test','Symbolic_Toolbox');
    if (res_partial(1)==0)
        disp('symbolic toolbox missing!');
    end

    % check if optimization toolbox is available
    res_partial(2) = license('test','Optimization_Toolbox');
    if (res_partial(2)==0)
        disp('optimization toolbox missing!');
    end

    p = path;

    % check if MPT toolbox is available
    res_partial(3) = contains(p,'mpt');
    if (res_partial(3)==0)
        disp('MPT toolbox missing!');
    end

    % check if YALMIP toolbox is available
    res_partial(4) = contains(p,'yalmip');
    if (res_partial(4)==0)
        disp('YALMIP toolbox missing!');
    end

    % check if CORA is on the path
    p = path;
    res_partial(5) = contains(p,'contDynamics'); 
    res_partial(6) = contains(p,'contSet'); 
    res_partial(7) = contains(p,'hybridDynamics');

    res = all(res_partial);

%------------- END OF CODE --------------
