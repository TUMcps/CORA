 %choose projection and plot------------------------------------------------
 options.projectedDimensions = [1 2];
 options.plotType = 'b';
 plot(PHA,'reachableSet',options); %plot reachable set
 plot(options.R0,options.projectedDimensions,'blackFrame'); %plot initial set
 plot(PHA,'simulation',options); %plot simulation
 %--------------------------------------------------------------------------