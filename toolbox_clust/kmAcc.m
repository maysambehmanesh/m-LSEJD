function [a,labels] = kmAcc(V,clusters,classes)

% note: eigenspaces can be first normalized to have rows with norm 1
V = bsxfun(@times, V, 1./(sum(V.^2, 2)));
labels=kmeans(V(:,1:clusters),clusters,'Replicates',100);
a=accuracy(labels,classes);

end