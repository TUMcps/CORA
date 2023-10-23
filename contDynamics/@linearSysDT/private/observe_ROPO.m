function [R,tcomp] = observe_ROPO(obj,options)
% observe_ROPO - computes the guaranteed state estimation approach
%    according to the set membership approach, see [1]. This is a strip
%    method.
%
% Syntax:
%    [R,tcomp] = observe_ROPO(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] Vicino, A., & Zappa, G. (1996). Sequential approximation of 
%        feasible parameter sets for identification with set membership 
%        uncertainty. IEEE Transactions on Automatic Control, 41(6), 774-785.
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Carlos Valero
% Written:       09-March-2021
% Last update:   09-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic
%time period
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% width of strips
options.sigma = supremum(abs(interval(options.V)));


% Intersection
y = options.y(:,1);
R{1} = options.R0;
Rnext=R{1};
% loop over all time steps
for k = 1:length(tVec)-1
    % Intersection
    y = options.y(:,k+1);
    % Update uncertain input
    Uadd = obj.B*options.uTransVec(:,k) + options.W;
    % Prediction
    Rnext = obj.A*Rnext+Uadd+obj.c;
    % Reduction to parallelotope
    Rnext=reduce(Rnext,options.reductionTechnique,1);
    % Strips Contruction
    p0=(obj.C./options.sigma)';
    c0=y./options.sigma;
    T = Rnext.generators;
    thetac = Rnext.center;
    no=size(p0,2);
    for i=1:no
        [T, thetac,~,~]=aux_Parallelotope(T,thetac,c0(i),p0(:,i));
    end
    Rnext=zonotope([thetac,T]);
    % Store result
    R{k+1} = Rnext;
end
tcomp=toc;
end


% Auxiliary functions -----------------------------------------------------

function [T, thetac,p0,c0]=aux_Parallelotope(t,thetac,c0,p0)
    %This function make one iteration and find the 
    %minimum parallelotope for n+1 strips
    %t E to R(nxn)
    %thetac E R(nx1)
    %c0 E R(1x1)
    %p0 E R(nx1)
    P0=0;           %initializing Sum P0=sum_{i=1}^{n}{p0*(ti)}
    n=size(t,1);    %Size of the parallelotope   
    %First it is necessary ensure that p0*(ti)>=0 for all i
    %if it doesn't change ti for -ti
    t0=p0/(norm(p0,2)^2);
    %%%% Now let's just for knowing who is E0p and E0n
    for i=1:n
        if p0'*t(:,i)<0
            t(:,i)=-t(:,i);
        end
        P0=P0+p0'*t(:,i);
    end
    E0p=(p0'*thetac-c0)+P0;     
    E0n=(p0'*thetac-c0)-P0;
%     if (E0p < 0) || (E0n > 0)
%         T=t;
%         return;
%     end
    %%%% First Stage. Tightening the measurement strip S0
    rop=min(1,E0p);
    ron=min(1,-E0n);
    p1=p0;
    %%%%%%%%% defining the new strip So
    p0=2*p0/(rop+ron);      
    c0=(2/(rop+ron))*(c0+0.5*(rop-ron));

    %%%% Second Stage. Reducing the parallelotope
    %%%% Let's start following the indications of Vicino 1996
    %Pre allocating Memory variables
    rp=zeros(1,n);
    rn=zeros(1,n);
    T=zeros(n);
    Imax=zeros(1,n+1);
    %%%% Building a new Parallelotope
    thetaC=0;
    for i=1:n
        K=p1'*t(:,i);
        if K==0
            rp(1,i)=1;
            rn(1,i)=1;
        else
            rp(1,i)=min(1,((1-E0n)/K)-1); 
            rn(1,i)=min(1,((1+E0p)/K)-1);
        end
        thetaC= thetaC+(rp(1,i)-rn(1,i))*t(:,i);
        t(:,i)=0.5*(rp(1,i)+rn(1,i))*t(:,i);
    end
    thetac=thetac+0.5*thetaC;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Third Stage. Selecting the minimal volume parallelotope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t=[t0 t];
    for i=1:n+1
        Imax(1,i)=p0'*t(:,i);
    end
    %our new t0 will be i=1 therefore
    [~,index]=max(Imax);
    if index==1
        T=t(:,2:n+1); 
    else
        k=1/(p0'*t(:,index));
        thetac=thetac+(k)*t(:,index)*(c0-(p0'*thetac));
        for i=2:n+1
            if i==index
                T(:,i-1)=k*t(:,index);
            else
                T(:,i-1)=t(:,i)-k*(p0'*t(:,i))*t(:,index);
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
