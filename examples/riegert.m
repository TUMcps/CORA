
close all; clear all;

load zonotope_test_2.mat
% load Polytope.mat

z = zonotope(rand(3,1),rand(3,5));

% Z = zonotope(zeros(3,1),rand(3,3));
% A = rand(6,3);
% Z2 = A*Z;
% G2 = generators(Z2);
% for i=1:5
% 	G2(:,end+1) = rand(1)*G2(:,1) + rand(1)*G2(:,2) + rand(1)*G2(:,3);
% end
% z = zonotope(center(Z2),G2);

% project to subspace, then convert
% G = generators(z);
% [U,~,~] = svd(G);
% G_3D = U'*G*U;
% G_3D = G_3D(1:3,:);
% Z_3D = zonotope(zeros(3,1),G_3D);
% P_3D = polytope(Z_3D);
% figure;
% projDim = {[1,2],[1,3],[2,3]};
% for i=1:length(projDim)
%     subplot(1,3,i); hold on; box on;
%     % plot sets
%     h_z = plot(Z_3D,projDim{i},'b');
%     h_p = plot(polytope(Z_3D),projDim{i},'r--');
%     % labels, legend
%     xlabel("x_" + projDim{i}(1));
%     ylabel("x_" + projDim{i}(2));
% end

% conversions
% ...using constructor of mptPolytope class
P_mptPolytope = mptPolytope(z);
% ...using @zonotope/polytope.m
P_polytope = polytope(z);
% ...using vertices of z and then instantiating mptPolytope object
V = vertices(z);
P_vertices = mptPolytope(V');

% visualization
projDim = {[1,2],[1,3],[2,3]};
% projDim = {[1,2],[3,4],[5,6]};
figure;
for i=1:length(projDim)
    subplot(1,3,i); hold on; box on;
    % plot sets
    h_z = plot(z,projDim{i},'b');
    h_mpt = plot(P_mptPolytope,projDim{i},'k.');
    h_p = plot(P_polytope,projDim{i},'g--');
    h_v = plot(P_vertices,projDim{i},'rx');
    % labels, legend
    xlabel("x_" + projDim{i}(1));
    ylabel("x_" + projDim{i}(2));
    legend([h_z,h_mpt,h_p,h_v],'Original zonotope','Conversion with mptPolytope',...
        'Conversion with polytope()','Conversion via vertices');
end