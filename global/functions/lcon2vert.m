function [V]=lcon2vert(A,b,Aeq,beq,TOL,checkbounds)
% lcon2vert - compute the vertices of a polytope
%
% Syntax:  
%    V = lcon2vert(A,b,Aeq,beq,TOL,checkbounds)
%
% Inputs:
%    A - inequality constraint matrix ( A*x < b )
%    b - inequality constraint vector ( A*x < b )
%    Aeq - equality constraint matrix ( A*x = b )
%    beq - equality constraint vector ( A*x = b )
%    TOL - tolerance for vertex computation
%    checkbounds - check if the polytope is bounded before the vertices are
%                  computed (0 or 1)
%
% Outputs:
%    V - matrix containing the vertices
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      09-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


  % initial argument parsing
  if nargin<5 || isempty(TOL), TOL=1e-10; end
  if nargin<6, checkbounds=true; end
  
  switch nargin 
      
      case 0          
            error('At least 1 input argument required');

      case 1
            b=[]; Aeq=[]; beq=[]; 
        
      case 2       
            Aeq=[]; beq=[];
          
      case 3      
            error('Since argument Aeq specified, beq must also be specified');        
  end
  
  b=b(:); beq=beq(:);
  
  if xor(isempty(A), isempty(b)) 
     error('Since argument A specified, b must also be specified');
  end
      
  if xor(isempty(Aeq), isempty(beq)) 
     error('Since argument Aeq specified, beq must also be specified');
  end
  
  nn=max(size(A,2)*~isempty(A),size(Aeq,2)*~isempty(Aeq));
  
  if ~isempty(A) && ~isempty(Aeq) && ( size(A,2)~=nn || size(Aeq,2)~=nn)    
     error('A and Aeq must have the same number of columns if both non-empty');  
  end
  
  
 % normalize constraints
 [A,b]=rownormalize(A,b);
 [Aeq,beq]=rownormalize(Aeq,beq);
 
 % extract all inequality constraints that describe equality constraints
 [A,b,Aeq,beq] = extractEqConstr(A,b,Aeq,beq,TOL);
 
 inequalityConstrained=~isempty(A);  
 equalityConstrained=~isempty(Aeq);

  
 % solve equality constraints        
 if equalityConstrained

   Neq=null(Aeq);   
   x0=pinv(Aeq)*beq;

   if norm(Aeq*x0-beq)>TOL*norm(beq)  % infeasible

      V=[]; 
      return;

   elseif isempty(Neq)

       if inequalityConstrained && ~all(A*x0<=b)

          V=[]; 
          return;              

       else % inequality constraints all satisfied, including vacuously

           V=x0(:).';           
           return
       end
   end
 end  
   
 % transform inequality constraints to the null-space from equality constr.
 if inequalityConstrained && equalityConstrained
     
   AAA=A*Neq;
   bbb=b-A*x0;
    
 elseif inequalityConstrained
      
    AAA=A;
    bbb=b;
   
 elseif equalityConstrained && ~inequalityConstrained
       error('Non-bounding constraints detected. (Consider box constraints on variables.)')     
 end
  
 % remove redundant constraints
 poly = Polyhedron(AAA,bbb);
 poly.minHRep();
 AAA = poly.A;
 bbb = poly.b;
  
  
 % calculat the vertices from the inequality constraints
 nnn=size(AAA,2);
  
 if nnn==1      % Special case
      
     idxu=sign(AAA)==1;
     idxl=sign(AAA)==-1;
     idx0=sign(AAA)==0;
     
     Q=bbb./AAA;
     U=Q; 
     U(~idxu)=inf;
     L=Q;
     L(~idxl)=-inf;

     [ub,uloc]=min(U);
     [lb,lloc]=max(L);
     
     if ~all(bbb(idx0)>=0) || ub<lb %infeasible
         
         V=[]; nr=[]; nre=[];
         return
         
     elseif ~isfinite(ub) || ~isfinite(lb)      
         error('Non-bounding constraints detected. (Consider box constraints on variables.)')
     end
      
     Zt=[lb;ub];
     
     if nargout>1
        nr=unique([lloc,uloc]); nr=nr(:);
     end   
      
 else           % Convert inequality constraints   

     Zt=con2vert(AAA,bbb,TOL,checkbounds); 
 end
  

 % convert vertices from null-space back to original space
  if equalityConstrained && ~isempty(Zt)     
      V=bsxfun(@plus,Zt*Neq.',x0(:).');      
  else     
      V=Zt;  
  end
  
  
  
  
 % AUXILIARY FUNCTIONS ----------------------------------------------------

 function [V,nr] = con2vert(A,b,TOL,checkbounds)
% CON2VERT - convert a convex set of constraint inequalities into the set
%            of vertices at the intersections of those inequalities;i.e.,
%            solve the "vertex enumeration" problem. Additionally,
%            identify redundant entries in the list of inequalities.
% 
% V = con2vert(A,b)
% [V,nr] = con2vert(A,b)
% 
% Converts the polytope (convex polygon, polyhedron, etc.) defined by the
% system of inequalities A*x <= b into a list of vertices V. Each ROW
% of V is a vertex. For n variables:
% A = m x n matrix, where m >= n (m constraints, n variables)
% b = m x 1 vector (m constraints)
% V = p x n matrix (p vertices, n variables)
% nr = list of the rows in A which are NOT redundant constraints
% 
% NOTES: (1) This program employs a primal-dual polytope method.
%        (2) In dimensions higher than 2, redundant vertices can
%            appear using this method. This program detects redundancies
%            at up to 6 digits of precision, then returns the
%            unique vertices.
%        (3) Non-bounding constraints give erroneous results; therefore,
%            the program detects non-bounding constraints and returns
%            an error. You may wish to implement large "box" constraints
%            on your variables if you need to induce bounding. For example,
%            if x is a person's height in feet, the box constraint
%            -1 <= x <= 1000 would be a reasonable choice to induce
%            boundedness, since no possible solution for x would be
%            prohibited by the bounding box.
%        (4) This program requires that the feasible region have some
%            finite extent in all dimensions. For example, the feasible
%            region cannot be a line segment in 2-D space, or a plane
%            in 3-D space.
%        (5) At least two dimensions are required.
%        (6) See companion function VERT2CON.
%        (7) ver 1.0: initial version, June 2005
%        (8) ver 1.1: enhanced redundancy checks, July 2005
%        (9) Written by Michael Kleder
%
%Modified by Matt Jacobson - March 30, 2011
% 


   % Check if the polytope is bounded
   if checkbounds
       
        [~,bb,~,bbeq]=vert2lcon(A,TOL);

        if any(bb<=0) || ~isempty(bbeq)
            error('Non-bounding constraints detected. (Consider box constraints on variables.)')
        end

        clear bb bbeq
   end
 
   
   
   % Initialization -------------------------------------------------------
   
   % determine a point that is located inside the polytope and that is far
   % enough away from the facets of the polytope
   
   dim=size(A,2);
   
   if strictinpoly(b,TOL)  
       
       c=zeros(dim,1);
   
   else
            
        % slackfun > 0 => point c is located inside the polytope
        slackfun=@(c)b-A*c;

        % Initializer 0
        c = pinv(A)*b; 
        s=slackfun(c);

        % Initializer 1
        if ~approxinpoly(s,TOL) 

            c=Initializer1(TOL,A,b,c);
            s=slackfun(c);
        end

        % Attempt refinement (Initializer 2)
        if  ~approxinpoly(s,TOL)  

            c=Initializer2(TOL,A,b,c);
            s=slackfun(c);
        end

        % Attempt refinement (Chebychef center)
        if ~approxinpoly(s,TOL) 
            
            poly = Polyhedron(A,b);
            temp = chebyCenter(poly);
            c = temp.x;
            s = slackfun(c);
        end

        % No point inside the polytope could be cound
        if ~approxinpoly(s,TOL)            
            error('Unable to locate a point near the interior of the feasible region.');
        end


        % Refinement: If the determined point is located too close to the
        % polyotpe surface, push it to the interior to increase the
        % numerical stability
        if ~strictinpoly(s,TOL) 

            % determine the surfaces that are too close to the point c
            idx=(  abs(s)<=max(s)*TOL );

            % Treat the close surfaces as equality constraints
            Amod=A; bmod=b; 
            Amod(idx,:)=[]; 
            bmod(idx)=[];

            Aeq=A(idx,:);
            beq=b(idx);

            % Calculate the vertices to the newly generated polytope which
            % is a subspace of the original polytope
            faceVertices=lcon2vert(Amod,bmod,Aeq,beq,TOL,1);
            
            if isempty(faceVertices)
               error('Something''s wrong. Couldn''t find face vertices. Possibly polyhedron is unbounded.');
            end

            % loop over all possible vertices (choose a feasible vertex)
            foundSolution = 0;

            for i = 1:size(faceVertices,1)
                
                try
                    
                    % find the local recession cone vector 
                    c=faceVertices(i,:).';
                    s=slackfun(c);

                    idx=(  abs(s)<=max(s)*TOL );

                    Asub=A(idx,:); bsub=b(idx,:);

                    [aa,bb,aaeq,bbeq]=vert2lcon(Asub);
                    aa=[aa;aaeq;-aaeq];
                    bb=[bb;bbeq;-bbeq];

                    clear aaeq bbeq

                    [bmin,idx]=min(bb);

                     if bmin>=-TOL
                       error('Something''s wrong. We should have found a recession vector (bb<0).');
                     end      


                    % find intersection of polytope with line through facet centroid.
                    Aeq2=null(aa(idx,:)).';
                    beq2=Aeq2*c;  

                    linetips = lcon2vert(A,b,Aeq2,beq2,TOL,1);

                    if size(linetips,1)<2
                       error('Failed to identify line segment through interior. Possibly {x: Aeq*x=beq} has weak intersection with interior({x: Ax<=b}).');
                    end

                    % Take midpoint to the line as the refined point
                    lineCentroid=mean(linetips);

                    clear aa bb

                    c=lineCentroid(:);
                    s=slackfun(c);

                    % suitable point was found => end search
                    foundSolution = 1;
                    break;
                end
            end

            if ~foundSolution
               error('Could not determine a point inside the polytope!'); 
            end
        end

        b = s;
   end

   
   % Calculate Vertices ---------------------------------------------------
   
   % Normalize the constraints (unit vector b)
   D=bsxfun(@rdivide,A,b); 
    
   % Determine the combinations of inequalies that form a vertex
   try
        k = convhulln(D);
   catch 
        k = convhulln(D,{'Qs'});
   end
   nr = unique(k(:));
    
   % For each vertex, combine the inequalities that form the vertex and
   % solve a system of linear equations to determine the point
    
   G  = zeros(size(k,1),dim);
   ee=ones(size(k,2),1);
   discard=false( 1, size(k,1) );
    
   for ix = 1:size(k,1)
        
        F = D(k(ix,:),:);
        if lindep(F,TOL)<dim
            discard(ix)=1;
            continue; 
        end

        G(ix,:)=F\ee;     
   end
    
   G(discard,:)=[];
   
   % shift the vertices by the initialization point (back to original
   % space)
   V = bsxfun(@plus, G, c.'); 
    
   % discard all vertices that are identical
   [~,I]=unique( round(V*1e6),'rows');
   V=V(I,:);
    
return


function [c,fval]=Initializer1(TOL, A,b,c,maxIter)
% Try to determine a point inside the polytope by iterative minimization of
% objective function

    thresh=-10*max(eps(b));
    
    if nargin>4
     [c,fval]=fminsearch(@(x) max([thresh;A*x-b]), c,optimset('MaxIter',maxIter));
    else
     [c,fval]=fminsearch(@(x) max([thresh;A*x-b]), c); 
    end
    
return          


function c=Initializer2(TOL,A,b,c)
% Try to determine a point inside the polytope by iterative minimization of
% objective function
    
    maxIter=10000;
    [mm,~]=size(A);
    
    Ap=pinv(A);        
    Aaug=speye(mm)-A*Ap;
    Aaugt=Aaug.';
    
    M=Aaugt*Aaug;
    C=sum(abs(M),2);
    C(C<=0)=min(C(C>0));
    
    slack=b-A*c;
    slack(slack<0)=0;
     
    IterThresh=maxIter; 
    s=slack; 
    ii=0;
    
    while ii<=2*maxIter 
        
       ii=ii+1; 
       if ii>IterThresh
           IterThresh=IterThresh+maxIter;
       end          
          
       s=s-Aaugt*(Aaug*(s-b))./C;   
       s(s<0)=0;
       c=Ap*(b-s);
    end
   
return 




function [r,idx,Xsub]=lindep(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [r,idx,Xsub]=lindep(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% r: rank estimate
% idx:  Indices (into X) of linearly independent columns
% Xsub: Extracted linearly independent columns of X

   if ~nnz(X) % X has no non-zeros and hence no independent columns
       
       Xsub=[]; idx=[];
       return
   end

   if nargin<2, tol=1e-10; end
   

           
     [~, R, E] = qr(X,0); 
     
     diagr = abs(diag(R));


     % Rank estimation
     r = find(diagr >= tol*diagr(1), 1, 'last'); % rank estimation

     if nargout>1
      idx=sort(E(1:r));
        idx=idx(:);
     end
     
     
     if nargout>2
      Xsub=X(:,idx);                      
     end                     

     
 function [A,b]=rownormalize(A,b)
 % Modifies A,b data pair so that norm of rows of A is either 0 or 1
 
  if isempty(A), return; end
 
  normsA=sqrt(sum(A.^2,2));
  idx=normsA>0;
  A(idx,:)=bsxfun(@rdivide,A(idx,:),normsA(idx));
  b(idx)=b(idx)./normsA(idx);       
        
 function tf=approxinpoly(s,TOL)
 % Determines if a point is approximateley (up to tolerance) located inside 
 % the polytope
     
   smax=max(s);
   
   if smax<=0
      tf=false; return 
   end
   
   tf=all(s>=-smax*TOL);
   
  function tf=strictinpoly(s,TOL)
      
   smax=max(s);
   
   if smax<=0
      tf=false; return 
   end
   
   tf=all(s>=smax*TOL);
   
function [A,b,Aeq,beq]=vert2lcon(V,tol)
%An extension of Michael Kleder's vert2con function, used for finding the 
%linear constraints defining a polyhedron in R^n given its vertices. This 
%wrapper extends the capabilities of vert2con to also handle cases where the 
%polyhedron is not solid in R^n, i.e., where the polyhedron is defined by 
%both equality and inequality constraints.
% 
%SYNTAX:
%
%  [A,b,Aeq,beq]=vert2lcon(V,TOL)
%
%The rows of the N x n matrix V are a series of N vertices of a polyhedron
%in R^n. TOL is a rank-estimation tolerance (Default = 1e-10).
%
%Any point x inside the polyhedron will/must satisfy
%  
%   A*x  <= b
%   Aeq*x = beq
%
%up to machine precision issues.
%
%
%EXAMPLE: 
%
%Consider V=eye(3) corresponding to the 3D region defined 
%by x+y+z=1, x>=0, y>=0, z>=0.
%
% 
%   >>[A,b,Aeq,beq]=vert2lcon(eye(3))
%
%
%     A =
% 
%         0.4082   -0.8165    0.4082
%         0.4082    0.4082   -0.8165
%        -0.8165    0.4082    0.4082
% 
% 
%     b =
% 
%         0.4082
%         0.4082
%         0.4082
% 
% 
%     Aeq =
% 
%         0.5774    0.5774    0.5774
% 
% 
%     beq =
% 
%         0.5774


  % initialization
  
  if nargin<2, tol=1e-10; end

  [M,N]=size(V);
    
  if M==1
      A=[];b=[];
      Aeq=eye(N); beq=V(:);
      return
  end

  p=V(1,:).';
  X=bsxfun(@minus,V.',p);
    
  % In the following, we need Q to be full column rank 
  % and we prefer E compact.
    
  if M>N  % X is wide
        
     [Q, R, E] = qr(X,0);  % economy-QR ensures that E is compact.
                           % Q automatically full column rank since X wide
                           
  else % X is tall, hence non-solid polytope
        
     [Q, R, P]=qr(X);  % non-economy-QR so that Q is full-column rank.
     
     [~,E]=max(P);  % No way to get E compact. This is the alternative. 
        clear P
  end 
    
  diagr = abs(diag(R));

    
  if nnz(diagr)    
        
        % Rank estimation
        r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation   
    
        iE=1:length(E);
        iE(E)=iE;
       
        Rsub=R(1:r,iE).';

        if r>1

          [A,b]=vert2con(Rsub,tol);
         
        elseif r==1
            
           A=[1;-1];
           b=[max(Rsub);-min(Rsub)];
        end

        A=A*Q(:,1:r).';
        b=bsxfun(@plus,b,A*p);
        
        if r<N
           Aeq=Q(:,r+1:end).';      
           beq=Aeq*p;
        else
           Aeq=[];
           beq=[];
        end

   else % Rank=0. All points are identical
      
       A=[]; b=[];
       Aeq=eye(N);
       beq=p;
   end
   
           
           
function [A,b] = vert2con(V,tol)
% VERT2CON - convert a set of points to the set of inequality constraints
%            which most tightly contain the points; i.e., create
%            constraints to bound the convex hull of the given points
%
% [A,b] = vert2con(V)
%
% V = a set of points, each ROW of which is one point
% A,b = a set of constraints such that A*x <= b defines
%       the region of space enclosing the convex hull of
%       the given points
%
% For n dimensions:
% V = p x n matrix (p vertices, n dimensions)
% A = m x n matrix (m constraints, n dimensions)
% b = m x 1 vector (m constraints)
%
% NOTES: (1) In higher dimensions, duplicate constraints can
%            appear. This program detects duplicates at up to 6
%            digits of precision, then returns the unique constraints.
%        (2) See companion function CON2VERT.
%        (3) ver 1.0: initial version, June 2005.
%        (4) ver 1.1: enhanced redundancy checks, July 2005
%        (5) Written by Michael Kleder, 
%
%Modified by Matt Jacobson - March 29,2011
% 

try
    k = convhulln(V);
catch
    k = convhulln(V,{'Qs'});
end
c = mean(V(unique(k),:));


V = bsxfun(@minus,V,c);
A  = nan(size(k,1),size(V,2));

dim=size(V,2);
ee=ones(size(k,2),1);
rc=0;

for ix = 1:size(k,1)
    F = V(k(ix,:),:);
    if lindep(F,tol) == dim
        rc=rc+1;
        A(rc,:)=F\ee;
    end
end

A=A(1:rc,:);
b=ones(size(A,1),1);
b=b+A*c';

% eliminate duplicate constraints:
[A,b]=rownormalize(A,b);
[~,I]=unique( round([A,b]*1e6),'rows');

A=A(I,:); % NOTE: rounding is NOT done for actual returned results
b=b(I);
return
         
         
function [A,b,Aeq,beq] = extractEqConstr(A,b,Aeq,beq,tol)
% Extract all inequality constraints that combined form an equality
% constraint and replace them by the respective equality constraint
    
    % Heuristic: sort the constraints according to a hash function to make
    % it easier to identify similar constraints
    A_ = abs(A);
    temp = 1:size(A,2);
    [A_,ind] = sortrows([A_*temp',A_]);
    A = A(ind,:);
    b = b(ind,:);
    
    
    % loop over all constraints
    i = 1;
    while i < size(A,1)
        
       % get index of last row that has the same hash value as the current 
       % row
       lInd = i;
       for k = i+1:size(A,1)
           if A_(k,1) ~= A_(i,1)
              break; 
           else
              lInd = lInd + 1;
           end
       end
       
       % compare all rows that have the same hash value to detect possible
       % equality constraints
       indRem = [];
       
       for j = i:lInd
          for k = j+1:lInd
              if all(abs(A(j,:)+A(k,:)) < tol) && abs(b(j)+b(k)) < tol
                 indRem = [indRem;j;k];
                 Aeq = [Aeq;A(j,:)];
                 beq = [beq;b(j)];
              end
          end
       end
       
       indRem = unique(indRem);
       A(indRem,:) = [];
       b(indRem) = [];
       A_(indRem,:) = [];
       
       i = lInd + 1 - length(indRem);
    end 
 
%------------- END OF CODE --------------