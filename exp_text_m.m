setLocalPaths

clear,clc

% load full coupling results (used as a possible initialization)
load('data/text_results.mat','VV','EE','V','E','clusters','classes','a','NMI');

% make sure classes is an nx1 vector
classes = classes(:); 

opt.numTests = 5;       % the result will be averaged over these tests
opt.K = 20;             % the number of eigenvectors to consider as input
opt.percentages = 40;   % the percentage of point whose correspondences are known
opt.maxIter = 100;      % max number of iterations for the solver
opt.lambda = 1e4;       % matching term
opt.gamma = 1e1;        % mismatching term
opt.saveResults = 0;	% save the results (provide the filename or resuts.mat by default)
%opt.saveFilename = sprintf('./results_CD_p%d_K%d',opt.percentages(1),opt.K);
opt.alpha=1;            % impact coefficients of additional matching samples
opt.mNumber=4;          % Number of Neighbor for Local PCA
opt.topEign=4;          % Number of top eigenvectors in Local PCA
opt.neigh=1;            % neighbor
opt.neighNumber=1;      % number of nieghbors

% you can specify a previous set of experiments to load the correspondences
% from (this allows to compare the performances of a different set of
% parameters with the same initialization data)
%opt.previousTests = load ('results.mat');

opt.algorithm = 'interior-point';
results = runSparseExperiments(opt, VV, EE, V, E, clusters, classes);
printMetrics(results);

% to compare with previous tests
%printMetrics(opt.previousTests.results);
