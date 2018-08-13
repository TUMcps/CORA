function [HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(dim)
% updated: 21-September-2011, MA

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
syms t

%combine symbolic variables to vectors
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27; x28; x29; x30; x31; x32; x33; x34; x35; x36; x37; x38; x39; x40; x41; x42; x43; x44; x45; x46; x47; x48; x49; x50; x51; x52; x53; x54; x55; x56; x57; x58; x59; x60; x61; x62; x63; x64; x65; x66; x67; x68; x69; x70; x71; x72; x73; x74; x75; x76; x77; x78; x79; x80; x81; x82; x83; x84; x85; x86; x87; x88; x89; x90; x91; x92; x93; x94; x95; x96; x97; x98; x99; x100; x101];                   
u = [u1; u2];

%compute jacobians for loc 1:
[f, c1] = fetchDynamics(dim,t,x,u,p);
A1 = double(jacobian(f, x(1:dim)));
B1 = double(jacobian(f, u));


%compute jacobians for loc 2:
p.kTmp = p.k;
p.k = 0;

[f, c2] = fetchDynamics(dim,t,x,u,p);
    
A2 = double(jacobian(f, x(1:dim)));
B2 = double(jacobian(f, u));
p.k = p.kTmp;

%compute jacobians for loc 3:
p.alphaTmp = p.alpha;
p.alpha = -p.alpha;

[f, c3] = fetchDynamics(dim,t,x,u,p);

A3 = double(jacobian(f, x(1:dim)));
B3 = double(jacobian(f, u));
p.alpha = p.alphaTmp;


%specify linear systems
linSys1 = linearSys('linearSys',A1,eye(dim));
linSys2 = linearSys('linearSys',A2,eye(dim));
linSys3 = linearSys('linearSys',A3,eye(dim));

%define distant
dist = 1e3;
eps = 1e-6;

%specify hyperplanes
%1st
n1 = [-1; zeros(dim-1,1)];
d1 = -p.alpha;
h1 = halfspace(n1,d1); 

n1b = [-1; zeros(dim-1,1)];
d1b = -p.alpha - 1e3;
h1b = halfspace(n1b,d1b); 

%2nd
n2 = [-1; zeros(dim-1,1)];
d2 = p.alpha;
h2 = halfspace(n2,d2); 

n2b = [-1; zeros(dim-1,1)];
d2b = p.alpha + 1e3;
h2b = halfspace(n2b,d2b); 

%loc1:
%invariant
inv = interval([p.alpha; -dist*ones(dim-1,1)], dist*ones(dim,1));
%guard sets
guard1 = h1b;
%resets
reset1.A = eye(dim);
reset1.b = zeros(dim,1);
%transitions
trans{1} = transition(guard1,reset1,2,'a','b'); %--> next loc: 2
%specify location
loc{1} = location('loc1',1,inv,trans,linSys1);

%loc2:
%invariant
inv = interval([-p.alpha - eps; -dist*ones(dim-1,1)], [p.alpha + eps; dist*ones(dim-1,1)]);
%guard set 1
guard1 = h1;
guard2 = h2b;
%reset 1
reset1.A=eye(dim);
reset1.b=zeros(dim,1);
%transition 1
trans=[];
trans{1}=transition(guard1,reset1,1,'a','b'); %--> next loc: 1
trans{2}=transition(guard2,reset1,3,'a','b'); %--> next loc: 3
%specify location
loc{2}=location('loc2',2,inv,trans,linSys2);

%loc3:
%invariant
inv = interval(-dist*ones(dim,1), [-p.alpha; dist*ones(dim-1,1)]);
%guard set 1
guard1 = h2;
%reset 1
reset1.A=eye(dim);
reset1.b=zeros(dim,1);
%transition 1
trans=[];
trans{1}=transition(guard1,reset1,2,'a','b'); %--> next loc: 2
%specify location
loc{3}=location('loc3',3,inv,trans,linSys3);

%specify hybrid automaton
HA=hybridAutomaton(loc);
end


function [f, c1] = fetchDynamics(dim,t,x,u,p)

    if dim==7
        f = powertrain7Eq(t,x,u,p);
        c1 = powertrain7Eq(0,zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==9
        f = powertrain9Eq(t,x,u,p);
        c1 = powertrain9Eq(0,zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==11
        f = powertrain11Eq(t,x,u,p);
        c1 = powertrain11Eq(0,zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==13
        f = powertrain13Eq(t,x,u,p);
        c1 = powertrain13Eq(0,zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==15
        f = powertrain15Eq(t,x,u,p);
        c1 = powertrain15Eq(0,zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==17
        f = powertrain17Eq(t,x,u,p);
        c1 = powertrain17Eq(0,zeros(length(x),1),zeros(length(u),1),p);
    elseif dim==21
        f = powertrain21Eq(t,x,u,p);
        c1 = powertrain21Eq(0,zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==31
        f = powertrain31Eq(t,x,u,p);
        c1 = powertrain31Eq(0,zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==41
        f = powertrain41Eq(t,x,u,p);
        c1 = powertrain41Eq(0,zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==51
        f = powertrain51Eq(t,x,u,p);
        c1 = powertrain51Eq(0,zeros(length(x),1),zeros(length(u),1),p);  
    elseif dim==61
        f = powertrain61Eq(t,x,u,p);
        c1 = powertrain61Eq(0,zeros(length(x),1),zeros(length(u),1),p); 
    elseif dim==101
        f = powertrain101Eq(t,x,u,p);
        c1 = powertrain101Eq(0,zeros(length(x),1),zeros(length(u),1),p);     
    end

end