clear;
close all;
load('centaur_hks.mat');
%%============read mesh from the off file============%%
fid=fopen('/Users/compume/Documents/MATLAB/COMPUTER_VISION_HW4/computer_vision_HW4/centaur1.off');
fgetl(fid);
nos = fscanf(fid, '%d %d  %d', [3 1]);
nopts = nos(1);
notrg = nos(2);
coord = fscanf(fid, '%g %g  %g', [3 nopts]);
coord = coord';
triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;
%%we have added 1 because the vertex indices start from 0 in vtk format
fclose(fid);
hold on;
plot3(coord(:,1),coord(:,2),coord(:,3),'g.');

hold off
descriptorgeodesic = zeros(nopts,nopts);
    [v,c] = size(triang);
    for i = 1:v
        descriptorgeodesic(triang(i,1),triang(i,2)) = sqrt(sum((coord(triang(i,2),:)-coord(triang(i,1),:)).^2));
        descriptorgeodesic(triang(i,1),triang(i,3)) = sqrt(sum((coord(triang(i,3),:)-coord(triang(i,1),:)).^2));
        descriptorgeodesic(triang(i,2),triang(i,3)) = sqrt(sum((coord(triang(i,3),:)-coord(triang(i,2),:)).^2));
        descriptorgeodesic(triang(i,2),triang(i,1)) = descriptorgeodesic(triang(i,1),triang(i,2));
        descriptorgeodesic(triang(i,3),triang(i,1)) = descriptorgeodesic(triang(i,1),triang(i,3));
        descriptorgeodesic(triang(i,3),triang(i,2)) = descriptorgeodesic(triang(i,2),triang(i,3));
    end
G = sparse(descriptorgeodesic);

%% program clustering
S = 1;
dist = zeros(10,nopts);
index = zeros(1,100);
[maxim,startind] = max(sum((centaur1_hks).^2,2));
index(1) = startind;
allindex = 1:3400;
dist = zeros(100,3400);
clusterindex = [startind];
dist(1,:) = graphshortestpath(G,index(1));
for i = 1:99
    J = setdiff(allindex,clusterindex);
    [maxim,indexmax] = max(sum(dist(1:i,J)));
    indexmax = J(indexmax);
    clusterindex(i+1) = indexmax;
    dist(i+1,:) = graphshortestpath(G,indexmax);
    [M,I] = min(dist(1:i+1,:),[],1);
    k = find(I == i+1);
    k = setdiff(k,clusterindex);
    [maxim,index1] = max(sum((centaur1_hks(k,:)).^2,2));
    index(i+1) = k(index1);
    dist(i+1,:) = graphshortestpath(G,index(i+1));
    clusterindex(i+1) = k(index1);
    
    
    
    
    
    
    
    
    
end
% colormap
[minimum,cluster] = min(dist,[],1);
cluster = clusterindex(cluster);
ofid = fopen('centaur1.vtk','w');
    fprintf(ofid, '# vtk DataFile Version 3.0\n');
    fprintf(ofid,'vtk output\n');
    fprintf(ofid,'ASCII\n');
    fprintf(ofid,'DATASET POLYDATA\n');
    fprintf(ofid,'POINTS %d float\n', nopts);
    fprintf(ofid,'%g %g %g\n', coord');
    fprintf(ofid,'POLYGONS %d %d\n', notrg, 4*notrg);
    fprintf(ofid,'3 %d %d %d\n', triang'-1);
    fprintf(ofid,'\n');
    fprintf(ofid,'POINT_DATA %d\n', nopts);
    fprintf(ofid,'SCALARS distance_from float\n');
    fprintf(ofid,'LOOKUP_TABLE default\n');
    fprintf(ofid,'%g\n',  cluster');
    fclose(ofid);
    
   


