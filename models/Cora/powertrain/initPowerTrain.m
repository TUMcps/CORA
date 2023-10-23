function [HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(dim)
% initPowerTrain - power train example described in Sec. 6 in [1]
%
% Syntax:  
%    [HA,Zcenter,Zdelta,B1,B2,B3,c1,c2,c3] = initPowerTrain(dim)
%
% Inputs:
%    dim - dimension of the model
%
% Outputs:
%    HA - hybridAutomaton object
%    Zcenter,Zdelta - center and generator for the initial zonotope
%    B1,B2,B3 - input matrices for the three locations
%    c1,c2,c3 - constant offsets for the three locations
%
% References:
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"

% Author:       Matthias Althoff
% Written:      21-September-2011
% Last update:  23-December-2019
% Last revision:---

%------------- BEGIN CODE --------------

%set parameters of the powertrain
[p,omega_ref_center,omega_ref_delta,x1_0,delta_x1_0,xy_0,delta_xy_0,delta_T_m0] = powertrainParameters();

% initial values
Zcenter = [x1_0;...
           p.T_m0;...
           0;...
           omega_ref_center;...
           0;...
           omega_ref_center;...
           p.i*omega_ref_center];
Zdelta =  [delta_x1_0;...
           delta_T_m0;...
           0;...
           omega_ref_delta; ...
           0; ...
           omega_ref_delta; ...
           p.i*omega_ref_delta];   
       
% add further dimensions
for iLoop = 1:((dim-7)/2)
    Zcenter(end+1:end+2,1) = [xy_0; ...
               omega_ref_center];
    Zdelta(end+1:end+2,1) =  [delta_xy_0; ...
               omega_ref_delta];
end


% %shrink initial zonotope
% Zinit = enlarge(Zinit,0.05*ones(dim,1));

%obtain linear system 
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 
syms x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31  
syms x32 x33 x34 x35 x36 x37 x38 x39 x40 x41
syms x42 x43 x44 x45 x46 x47 x48 x49 x50 x51
syms x52 x53 x54 x55 x56 x57 x58 x59 x60 x61
syms x62 x63 x64 x65 x66 x67 x68 x69 x70 x71
syms x72 x73 x74 x75 x76 x77 x78 x79 x80 x81
syms x82 x83 x84 x85 x86 x87 x88 x89 x90 x91
syms x92 x93 x94 x95 x96 x97 x98 x99 x100 x101
syms u1 u2 

%combine symbolic variables to vectors
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27; x28; x29; x30; x31; x32; x33; x34; x35; x36; x37; x38; x39; x40; x41; x42; x43; x44; x45; x46; x47; x48; x49; x50; x51; x52; x53; x54; x55; x56; x57; x58; x59; x60; x61; x62; x63; x64; x65; x66; x67; x68; x69; x70; x71; x72; x73; x74; x75; x76; x77; x78; x79; x80; x81; x82; x83; x84; x85; x86; x87; x88; x89; x90; x91; x92; x93; x94; x95; x96; x97; x98; x99; x100; x101];                   
u = [u1; u2];

%compute jacobians for loc 1:
[f, c1] = fetchDynamics(dim,x,u,p);
A1 = double(jacobian(f, x(1:dim)));
B1 = double(jacobian(f, u));


%compute jacobians for loc 2:
p.kTmp = p.k;
p.k = 0;

[f, c2] = fetchDynamics(dim,x,u,p);
    
A2 = double(jacobian(f, x(1:dim)));
B2 = double(jacobian(f, u));
p.k = p.kTmp;

%compute jacobians for loc 3:
p.alphaTmp = p.alpha;
p.alpha = -p.alpha;

[f, c3] = fetchDynamics(dim,x,u,p);

A3 = double(jacobian(f, x(1:dim)));
B3 = double(jacobian(f, u));
p.alpha = p.alphaTmp;


%specify linear systems
linSys1 = linearSys('linearSys1',A1,eye(dim));
linSys2 = linearSys('linearSys2',A2,eye(dim));
linSys3 = linearSys('linearSys3',A3,eye(dim));

%define distant
dist = 1e3;
eps = 1e-6;

%specify hyperplanes
n = [1;zeros(dim-1,1)];

%1st
h1 = conHyperplane(n,p.alpha); 

%2nd
h2 = conHyperplane(n,-p.alpha); 


%loc1:
%invariant
inv = polytope(-n',-p.alpha);
%guard sets
guard1 = h1;
%resets
reset1.A = eye(dim);
reset1.c = zeros(dim,1);
%transitions
trans = transition(guard1,reset1,2); %--> next loc: 2
%specify location
loc = location('loc1',inv,trans,linSys1);

%loc2:
%invariant
inv = polytope([n';-n'],[p.alpha;p.alpha]);
%guard set 1
guard1 = h1;
guard2 = h2;
%reset 1
reset1.A=eye(dim);
reset1.c=zeros(dim,1);
%transition 1
clear trans
trans(1)=transition(guard1,reset1,1); %--> next loc: 1
trans(2)=transition(guard2,reset1,3); %--> next loc: 3
%specify location
loc(2)=location('loc2',inv,trans,linSys2);

%loc3:
%invariant
inv = polytope(n',-p.alpha);
%guard set 1
guard1 = h2;
%reset 1
reset1.A=eye(dim);
reset1.c=zeros(dim,1);
%transition 1
clear trans
trans=transition(guard1,reset1,2); %--> next loc: 2
%specify location
loc(3)=location('loc3',inv,trans,linSys3);

%specify hybrid automaton
HA=hybridAutomaton(loc);
end


function [f, c1] = fetchDynamics(dim,x,u,p)

    if dim==7
        f = powertrain7Eq(x,u,p);
        c1 = powertrain7Eq(zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==9
        f = powertrain9Eq(x,u,p);
        c1 = powertrain9Eq(zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==11
        f = powertrain11Eq(x,u,p);
        c1 = powertrain11Eq(zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==13
        f = powertrain13Eq(x,u,p);
        c1 = powertrain13Eq(zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==15
        f = powertrain15Eq(x,u,p);
        c1 = powertrain15Eq(zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==17
        f = powertrain17Eq(x,u,p);
        c1 = powertrain17Eq(zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==21
        f = powertrain21Eq(x,u,p);
        c1 = powertrain21Eq(zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==31
        f = powertrain31Eq(x,u,p);
        c1 = powertrain31Eq(zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==41
        f = powertrain41Eq(x,u,p);
        c1 = powertrain41Eq(zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==51
        f = powertrain51Eq(x,u,p);
        c1 = powertrain51Eq(zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==61
        f = powertrain61Eq(x,u,p);
        c1 = powertrain61Eq(zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==101
        f = powertrain101Eq(x,u,p);
        c1 = powertrain101Eq(zeros(length(x),1),zeros(length(u),1),p);     
    end

end

%------------- END OF CODE --------------