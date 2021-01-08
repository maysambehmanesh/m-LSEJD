function [neighborhood,distance]=findNeighbor(X,m)

%input:   X: data  ,  K:number of neigbor
%output   index: closest neigbor    distance:

[distance,index]=sort(X);

neighborhood = index(2:(1+m),:);
distance=distance(2:(1+m),:);




end