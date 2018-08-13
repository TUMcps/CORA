function [matZ_A,matZ_B] = initRLCnonStiff()
% updated: 04-October-2010, MA


%nr of nodes
nrOfNodes=20;

%build matrix zonotope
A0=zeros(2*nrOfNodes);


%p_1
%p{1}=1/L;

%first node
tmpA=A0;
tmpA(nrOfNodes+1,1)=1;

zonB{1}=zeros(2*nrOfNodes,1);
zonB{1}(nrOfNodes+1,1)=-1;

%other nodes
for i=2:(nrOfNodes)
    tmpA(nrOfNodes+i,i-1)=-1;
    tmpA(nrOfNodes+i,i)=+1;
end
zonA{1}=tmpA;

%p_2
%p{2}=1/C;
tmpA=A0;
for i=1:(nrOfNodes-1)
    tmpA(i,nrOfNodes+i)=-1;
    tmpA(i,nrOfNodes+i+1)=+1;
end
tmpA(nrOfNodes,2*nrOfNodes)=-1;
zonA{2}=tmpA;

%p_3
%p{3}=Rdriver/L;
tmpA=A0;
tmpA(nrOfNodes+1,nrOfNodes+1)=-1;
zonA{3}=tmpA;

%p_4
%p{4}=R/L;
tmpA=A0;
for i=2:(nrOfNodes)
    tmpA(nrOfNodes+i,nrOfNodes+i)=-1;
end
zonA{4}=tmpA;

%good parameters
R = interval(0.99,1.01); 
Rdriver = interval(9.9,10.1);
L = interval(1e-10,1e-10); 
C = interval(3.99e-13,4.01e-13); 

%obtain parameter ranges
pVal(1)=1/L;
pVal(2)=1/C;
pVal(3)=Rdriver/L;
pVal(4)=R/L;

%get centers and normalized deltas
pCenter = mid(pVal);
pRad = rad(pVal);

%normalize cell matrices
%system matrix
Acenter = 0*zonA{1};
for i=1:length(zonA)
    %add center matrices
    Acenter = Acenter + pCenter(i)*zonA{i};
    Arad{i} = pRad(i)*zonA{i};
end

Bcenter = pCenter(1)*zonB{1};
Brad{1} = pRad(1)*zonB{1};

%instantiate matrix zonotopes
matZ_A=matZonotope(Acenter,Arad);
matZ_B=matZonotope(Bcenter,Brad);

%correction due to time scaling
matZ_A=1e-9*matZ_A;
matZ_B=1e-9*matZ_B;