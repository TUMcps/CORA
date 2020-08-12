function EstSet = SMA_SISO_Zonotope(Estimator,PpD,xh0,options,OffGain)

% To compute the SISO technique based state estimation USING SMA approach

ts = options.stepSize;  % Step size of simulation
tspan = 0:ts:options.T;   % total time span for run
%---------Initialization-----------%
EstSet.ZXh = xh0.ZXin; EstSet.ZXp = xh0.ZXin; EstSet.Xhp(:,1) = center(xh0.ZXin);EstSet.Xpp(:,1) = center(xh0.ZXin);
EstSet.XhH{1} = generators(xh0.ZXin);EstSet.XpH{1} = generators(xh0.ZXin);  
isWtCombastel =  (isfield(options,'RedTech') && strcmp(options.RedTech,'Wtcombastel'));
isVolMinII  =  (isfield(Estimator,'Name') && strcmp(Estimator.Name,'VolMinII'));
%-----------Start------%
%-----------Using generic Reduction Technique --------------%  
if  ~isWtCombastel
for k = 1:length(tspan)-2         % The coupled SMA based approach 
EstSet.ZXh(:,k) = reduce(zonotope([EstSet.Xhp(:,k),EstSet.XhH{k}]),options.RedTech,options.ZOrder);
EstSet.XhH{k} = generators(EstSet.ZXh(:,k));
EstSet.Xpp(:,k+1) = PpD.A*EstSet.Xhp(:,k) + PpD.B*PpD.U(:,k);% Center of predicted zonotope at k+1 instant
EstSet.XpH{k+1} =[PpD.A*EstSet.XhH{k}, PpD.E*generators(xh0.ZW)]; % Generators of predicted zonotope at time k+1 
EstSet.ZXp(:,k+1) = zonotope([EstSet.Xpp(:,k+1),EstSet.XpH{k+1}]); % The predicted zonotope
%-------------------------------------------------
for j = 1: size(PpD.C,1)
     if  j == 1
%-------------------         
if ~isVolMinII         
EstSet.OGain{k}(:,j) =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XpH{k+1},EstSet.Xpp(:,k+1),PpD.Yp(j,k+1),j,OffGain);   
EstSet.EstX{j} =  EstSet.Xpp(:,k+1) + EstSet.OGain{k}(:,j)*(PpD.Yp(j,k+1)-(PpD.C(j,:)*EstSet.Xpp(:,k+1)));     
EstSet.XhHt{j} = [(eye(size(PpD.A,1))- (EstSet.OGain{k}(:,j)*PpD.C(j,:)))*EstSet.XpH{k+1}, EstSet.OGain{k}(:,j)*PpD.F(j,:)*generators(xh0.ZV)]; 
elseif isVolMinII
EstSet.ZXhHt{j} =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XpH{k+1},EstSet.Xpp(:,k+1),PpD.Yp(j,k+1),j,OffGain);  
end
%------------------
    elseif j~=1
%-------------------         
if ~isVolMinII
EstSet.OGain{k}(:,j) =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XhHt{j-1},EstSet.EstX{j-1},PpD.Yp(j,k+1),j,OffGain); 
EstSet.EstX{j} =  EstSet.EstX{j-1} + EstSet.OGain{k}(:,j)*(PpD.Yp(j,k+1)-(PpD.C(j,:)*EstSet.EstX{j-1}));   
EstSet.XhHt{j} = [(eye(size(PpD.A,1))- (EstSet.OGain{k}(:,j)*PpD.C(j,:)))*EstSet.XhHt{j-1}, EstSet.OGain{k}(:,j)*PpD.F(j,:)*generators(xh0.ZV)];
%-----------------------
elseif isVolMinII
EstSet.ZXhHt{j} =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),generators(EstSet.ZXhHt{j-1}),center(EstSet.ZXhHt{j-1}),PpD.Yp(j,k+1),j,OffGain);      
end
     end
end 
%----------------------------------
if ~isVolMinII
EstSet.Xhp(:,k+1) = EstSet.EstX{end}  ;EstSet.XhH{k+1} =  EstSet.XhHt{end};  
elseif isVolMinII
EstSet.Xhp(:,k+1) = center(EstSet.ZXhHt{end}); EstSet.XhH{k+1} = generators(EstSet.ZXhHt{end});
end
%--------------------------
EstSet.ZXh(:,k+1) = zonotope([EstSet.Xhp(:,k+1),EstSet.XhH{k+1}]); 
end
end

%---------------Using Wt. Combastel reduction Operator
if  isWtCombastel
for k = 1:length(tspan)-2         % The coupled SMA based approach 
EstSet.ZXh(:,k) = WeightedCombastel(zonotope([EstSet.Xhp(:,k),EstSet.XhH{k}]),options.ZOrder,OffGain.IOA_W_NomG);
EstSet.XhH{k} = generators(EstSet.ZXh(:,k));
EstSet.Xpp(:,k+1) = PpD.A*EstSet.Xhp(:,k) + PpD.B*PpD.U(:,k);% Center of predicted zonotope at k+1 instant
EstSet.XpH{k+1} =[PpD.A*EstSet.XhH{k}, PpD.E*generators(xh0.ZW)]; % Generators of predicted zonotope at time k+1 
EstSet.ZXp(:,k+1) = zonotope([EstSet.Xpp(:,k+1),EstSet.XpH{k+1}]); % The predicted zonotope
%-------------------------------------------------
for j = 1: size(PpD.C,1)
     if  j == 1
%-------------------         
if ~isVolMinII         
EstSet.OGain{k}(:,j) =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XpH{k+1},EstSet.Xpp(:,k+1),PpD.Yp(j,k+1),j,OffGain);   
EstSet.EstX{j} =  EstSet.Xpp(:,k+1) + EstSet.OGain{k}(:,j)*(PpD.Yp(j,k+1)-(PpD.C(j,:)*EstSet.Xpp(:,k+1)));     
EstSet.XhHt{j} = [(eye(size(PpD.A,1))- (EstSet.OGain{k}(:,j)*PpD.C(j,:)))*EstSet.XpH{k+1}, EstSet.OGain{k}(:,j)*PpD.F(j,:)*generators(xh0.ZV)]; 
elseif isVolMinII
EstSet.ZXhHt{j} = Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XpH{k+1},EstSet.Xpp(:,k+1),PpD.Yp(j,k+1),j,OffGain);   
end
%------------------
    elseif j~=1
%-------------------         
if ~isVolMinII
EstSet.OGain{k}(:,j) =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),EstSet.XhHt{j-1},EstSet.EstX{j-1},PpD.Yp(j,k+1),j,OffGain); 
EstSet.EstX{j} =  EstSet.EstX{j-1} + EstSet.OGain{k}(:,j)*(PpD.Yp(j,k+1)-(PpD.C(j,:)*EstSet.EstX{j-1}));   
EstSet.XhHt{j} = [(eye(size(PpD.A,1))- (EstSet.OGain{k}(:,j)*PpD.C(j,:)))*EstSet.XhHt{j-1}, EstSet.OGain{k}(:,j)*PpD.F(j,:)*generators(xh0.ZV)];
%-----------------------
elseif isVolMinII
EstSet.ZXhHt{j} =  Estimator.ComputeGainSISO(PpD,PpD.C(j,:),PpD.F(j,:),generators(EstSet.ZXhHt{j-1}),center(EstSet.ZXhHt{j-1}),PpD.Yp(j,k+1),j,OffGain);      
end
     end
end 
%----------------------------------
if ~isVolMinII
EstSet.Xhp(:,k+1) = EstSet.EstX{end}  ;EstSet.XhH{k+1} =  EstSet.XhHt{end};  
elseif isVolMinII
EstSet.Xhp(:,k+1) = center(EstSet.ZXhHt{end}); EstSet.XhH{k+1} = generators(EstSet.ZXhHt{end});
end
%--------------------------
EstSet.ZXh(:,k+1) = zonotope([EstSet.Xhp(:,k+1),EstSet.XhH{k+1}]); 
end
end



end