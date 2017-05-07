%% CENTAUR1 K-MEAN, 
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


[maxim,startind] = max(sum((centaur1_hks).^2,2));
index = [startind];
allindex = 1:3400;
dist = graphallshortestpaths(G);

for i = 1:9
   if i == 1
       [maxim,k] = max(dist(index,:));
   else
       J = setdiff(allindex,index);
       [maxim,k] = max(sum(dist(index,J).^2,1));
       k = J(k);
       
   end
      
   index(i+1) = k ; 
   %[maxim1,cluster] = min(dist(index,:),[],1);
   %cluster = index(cluster);
   %k1 = find(cluster == index(i+1));
   %[ind,l] = max(sum((centaur1_hks(k1,:)).^2,2));
   %index(i+1) = k1(l);
end
% colormap
%[minimum,cluster] = min(dist(index,:),[],1);
%cluster = index(cluster);
%cluster(index) = 3400;
k = 1;
lamda = 0.8;

R1 =repmat(sqrt(sum(centaur1_hks.^2,2)),1,101);
R2 =repmat(sqrt(sum(coord.^2,2)),1,3);
centaur1_hks = centaur1_hks./R1;
coord1 = coord./R2;
measure1 = [lamda*centaur1_hks (1-lamda)*coord1];
cost1 = mean(measure1,2);
dist1 = zeros(10,3400);
H = 1;
for j = 1:500
    index1 = index;
    dist2 = dist(index,:);
  
  [minim,cluster] = min(dist2,[],1);
  minind = zeros(1,10);
  for i = 1:10
      k = find(cluster == i);
      L = length(k);
      meancoord = mean(measure1(k,:));
      [minim,minind(i)]   = min(sum((measure1(k,:) - repmat(meancoord,L,1)).^2,2));
      minind(i) = k(minind(i));
      index(i) = minind(i);
  end
  if (sum(index - index1) == 0)
      H = 0;
  end
  
    
end
colorvec = zeros(1,10);
for i = 1:10
    colorvec(i) =  200*(i-1);
end
cluster = colorvec(cluster);

ofid = fopen('centaur2.vtk','w');
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
    