function HA = gearbox_old(~)


%% Generated on 20-Apr-2018

%----------------Automaton created from Component 'mesh'-------------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (mesh.Clock_1 X mesh.Stateflow_2):
%  state x := [t; vx; vy; px; py; II]
%  input u := [uDummy]

%---------------Component mesh.Clock_1 X mesh.Stateflow_2------------------

%------------------------State loc01 X move_free---------------------------

%% equation:
%   t' == 1
%   &&
%   vx'==Fs/ms &
%   vy'==-Rs*Tf/Jg2 &
%   px'==vx &
%   py'==vy &
%   I'==0
dynA = ...
[0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,0,0,0];
dynB = ...
[0;0;0;0;0;0];
dync = ...
[1;21.875;-0.1142857143;0;0;0];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   t<=0.5
%   &&
%   px<=deltap &
%   py<=-px*0.726542528005361 & py>=px*0.726542528005361
invA = ...
[1,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0.726542528,1,0;0,0,0,0.726542528,-1,0];
invb = ...
[0.5;-0.003;-0;-0];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
%% equation:
%   I:=I+(vx*0.587785252292473+vy*0.809016994374947)*(zeta+1)*ms*mg2/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473)) &
%   vx:=(vx*(ms*0.809016994374947*0.809016994374947-mg2*zeta*0.587785252292473*0.587785252292473)+vy*(-(zeta+1)*mg2*0.587785252292473*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473)) &
%   vy:=(vx*(-(zeta+1)*ms*0.587785252292473*0.809016994374947)+vy*(mg2*0.587785252292473*0.587785252292473-ms*zeta*0.809016994374947*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))
resetA = ...
[1,0,0,0,0,0;0,-0.4232994906,-1.959003686,0,0,0;0,-0.3463431932,...
0.5232994906,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,7.748677518,10.66513964,0,0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   py>=-px*0.726542528005361 && 
%   vx*0.587785252292473+vy*0.809016994374947>0
guardA = ...
[0,0,0,-0.726542528,-1,0;0,-0.5877852523,-0.8090169944,0,0,0];
guardb = ...
[-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{1} = transition(guard, reset, 1, 'dummy', 'names');

%% equation:
%   I:=I+(vx*0.587785252292473-vy*0.809016994374947)*(zeta+1)*ms*mg2/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473)) &
%   vx:=(vx*(ms*0.809016994374947*0.809016994374947-mg2*zeta*0.587785252292473*0.587785252292473)+vy*((zeta+1)*mg2*0.587785252292473*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473)) &
%   vy:=(vx*((zeta+1)*ms*0.587785252292473*0.809016994374947)+vy*(mg2*0.587785252292473*0.587785252292473-ms*zeta*0.809016994374947*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))
resetA = ...
[1,0,0,0,0,0;0,-0.4232994906,1.959003686,0,0,0;0,0.3463431932,...
0.5232994906,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,7.748677518,-10.66513964,0,...
0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   py<=px*0.726542528005361 && 
%   vx*0.587785252292473-vy*0.809016994374947>0
guardA = ...
[0,0,0,-0.726542528,1,0;0,-0.5877852523,0.8090169944,0,0,0];
guardb = ...
[-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{2} = transition(guard, reset, 1, 'dummy', 'names');

%% equation:
%   I:=I+ms*vx+ms*vy &
%   vx:=0 &
%   vy:=0
resetA = ...
[1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,3.2,3.2,...
0,0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   px>=deltap & vx>=0 & vy>=0
guardA = ...
[0,0,0,-1,0,0;0,-1,0,0,0,0;0,0,-1,0,0,0];
guardb = ...
[0.003;-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{3} = transition(guard, reset, 2, 'dummy', 'names');

%% equation:
%   I:=I+ms*vx-ms*vy &
%   vx:=0 &
%   vy:=0
resetA = ...
[1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,3.2,-3.2,...
0,0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   px>=deltap & vx>=0 & vy<=0
guardA = ...
[0,0,0,-1,0,0;0,-1,0,0,0,0;0,0,1,0,0,0];
guardb = ...
[0.003;-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{4} = transition(guard, reset, 2, 'dummy', 'names');

%% equation:
%   I:=I-ms*vx+ms*vy &
%   vx:=0 &
%   vy:=0
resetA = ...
[1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,-3.2,3.2,...
0,0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   px>=deltap & vx<=0 & vy>=0
guardA = ...
[0,0,0,-1,0,0;0,1,0,0,0,0;0,0,-1,0,0,0];
guardb = ...
[0.003;-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{5} = transition(guard, reset, 2, 'dummy', 'names');

%% equation:
%   I:=I-ms*vx-ms*vy &
%   vx:=0 &
%   vy:=0
resetA = ...
[1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,-3.2,...
-3.2,0,0,1];
resetb = ...
[0;0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   px>=deltap & vx<=0 & vy<=0
guardA = ...
[0,0,0,-1,0,0;0,1,0,0,0,0;0,0,1,0,0,0];
guardb = ...
[0.003;-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{6} = transition(guard, reset, 2, 'dummy', 'names');

loc{1} = location('S1',1, inv, trans, dynamics);



%-------------------------State loc01 X meshed-----------------------------

%% equation:
%   t' == 1
%   &&
%   false
dynA = ...
[0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0];
dynB = ...
[0;0;0;0;0;0];
dync = ...
[1;0;0;0;0;0];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   t<=0.5
%   &&
%   
invA = ...
[1,0,0,0,0,0];
invb = ...
[0.5];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
loc{2} = location('S2',2, inv, trans, dynamics);



HA = hybridAutomaton(loc);


end