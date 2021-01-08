function A = adjacency_weighted_cosine(DATA, PARAM, t)

% Compute the weighted adjacency graph of the data set DATA with
% self-tuning gaussian kernel parameter. See 
% http://webee.technion.ac.il/~lihi/Demos/SelfTuningClustering.html
%
% A = adjacency(DATA, PARAM, t);
% 
% DATA - NxK matrix. Data points are rows. 
% PARAM - integer: number of nearest neighbors
% t - (typically) integer: parameter for gaussian kernel
%                          if not provided, self-tuning is applied
%
% Returns: A, sparse symmetric NxN matrix of distances between the
% adjacent points. 
%
% Example: 
% 
% A = adjacency(X,6) 
%   A contains the adjacency matrix for the data
%   set X. For each point, the distances to 6 adjacent points are
%   stored. Self-tuning is applied for gaussian kernel
%
% Note: the adjacency relation is symmetrized, i.e. if
% point a is adjacent to point b, then point b is also considered to be
% adjacent to point a.
%
%
% Original Author: 
%
% Mikhail Belkin 
% misha@math.uchicago.edu
%
% Modified By:
% Davide Eynard
% davide.eynard@usi.ch


if (nargin < 2) || ~isreal(PARAM)
  fprintf('ERROR: Too few arguments given or incorrect arguments.\n');
  return;
end

if (nargin < 3)
  t = 0;
end



n = size(DATA,1);
fprintf ('DATA: %d points in %d dimensional space.\n',n,size (DATA,2));
fprintf('Creating the adjacency matrix. Nearest neighbors, N=%d.\n', PARAM); 

  
A = sparse(n,n);

% NOTE: L2_distance is WAY FASTER than pdist!!! But here we need cosine...
%dt = L2_distance(DATA',DATA');
dt = squareform(pdist(DATA,'cosine'));
[Z,I] = sort (dt,2);

for i=1:n
    if ( mod(i,500) == 0) 
        fprintf('%d points processed.\n', i);
    end
    for j=2:PARAM+1
        if (t)
            % use gaussian kernel with t as a parameter
            A(i,I(i,j))= exp(-Z(i,j).^2/t);
        else
            % use self-tuning
            A(i,I(i,j))= exp(-Z(i,j).^2/(1e-8 + Z(i,PARAM+1)*Z(I(i,j),PARAM+1)));
            if isnan(A(i,I(i,j)))  
                % isnan might happen if neighbor nn+1 (plus all the
                % previous) has distance==0. To fix this I added the
                % 1e-8 above
                fprintf('[w] isNaN! i=%d,j=%d: %f, %f, %f\n',i,j,Z(i,j),Z(i,PARAM+1),Z(j,PARAM+1));
            end
        end
        A(I(i,j),i)= A(i,I(i,j));
    end    
end
    
end