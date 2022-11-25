function [Rnext,options] = post(obj,R,options)
% post - computes the reachable continuous probabilistic set for one time 
% step of a linear system
%
% Syntax:  
%    [Rnext,options] = post(obj,R,options)
%
% Inputs:
%    obj - linProbSys object
%    R - reachable set of the previous time step
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-September-2007 
% Last update:  08-September-2009
%               25-July-2016 (intervalhull replaced by interval)
%               16-July-2020 (update to integrate in CORA 2020)
% Last revision:---

%------------- BEGIN CODE --------------

%load homogeneous solution from options struct
Rhom=options.Rhom;
Rpar_det=options.Rpar_det;
Rpar_prob=options.Rpar_prob;
Raux_det=options.Raux_det;

%update 
eAt=obj.taylor.eAt;
pRinput=obj.taylor.pRinput;
Rhom=eAt*Rhom;
Rpar_det=Rpar_det+interval(Raux_det);
Rpar_prob=eAt*Rpar_prob+pRinput;
Raux_det=eAt*Raux_det;

%save homogeneous and particulate solution to options struct
options.Rhom=Rhom;
options.Rpar_det=Rpar_det;
options.Rpar_prob=Rpar_prob;
options.Raux_det=Raux_det;

%write results to reachable set struct Rnext
Rnext.tp=[];
Rnext.ti=Rhom+Rpar_prob+zonotope(Rpar_det);

%------------- END OF CODE --------------