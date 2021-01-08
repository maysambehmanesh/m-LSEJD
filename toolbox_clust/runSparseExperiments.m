function [results] = runSparseExperiments(opt, VV, EE, V, E, clusters, classes)

%- Default values

if ~exist( 'classes', 'var' ) || isempty(classes)
    fprintf('You need to provide classes (ground truth) for evaluation!\n\n');
    return;
end

if ~exist( 'clusters', 'var' ) || isempty(clusters)
    clusters = length(unique(classes));
    fprintf('[i] I automatically detected the number of clusters from the classes variable: %d\n',clusters);
end

if ~exist( 'V', 'var' ) || isempty(V)
    initJADE = 0;
    fprintf('[i] Variable V has not been provided: initialization is random\n');
else
    initJADE = 1;
    fprintf('[i] Variable V has been provided: initialization is JADE\n');
end

if ~isfield(opt, 'K')
    K = 20;
else
    K = opt.K;
end

if ~isfield(opt, 'Z')
    Z = clusters;
    fprintf('[i] Variable Z has not been provided: it has been set to #clusters=%d\n',clusters);
else
    Z = opt.Z;
end

if ~isfield(opt, 'numTests')
    opt.numTests = 10;
    fprintf('[i] Variable opt.numTests has not been provided: running %d tests\n',opt.numTests);
end

if ~isfield(opt, 'percentages')
    opt.percentages = [10 20 60 100];
end

if ~isfield(opt, 'maxIter')
    opt.maxIter = 100;
end

if ~isfield(opt, 'lambda')
    opt.lambda = 1e2;
end

if ~isfield(opt, 'gamma')
    opt.gamma = 1e0;
end

if ~isfield(opt, 'algorithm')
    opt.algorithm = 'sqp';
end

if ~isfield(opt, 'saveResults')
    opt.saveResults = 0;
end

if ~isfield(opt, 'saveFilename')
    opt.saveFilename = 'results';
end

if ~isfield(opt, 'previousTests')
    LOADCORRESP = 0;
    V = VV{1}; % used to calculate diffusion distances
    E = EE{1};
else
    LOADCORRESP = 1;
    % TODO: verify that the previous tests file is compatible
    % 1) it has to have the same number of percentages
    % 2) it has to have the same number of tests per percentage
    % (currently I just make sure I am passing the correct data files)
end

% Parameter
m=opt.mNumber;          
te=opt.topEign;          
s=opt.neigh;           
l=opt.neighNumber;           


% #points in modality 1 and 2
N = size( VV{1}, 1 );
M = size( VV{2}, 1 );

clear D
% getting the eigenvalues
for i = 1:length(VV)
    D(:,i) = diag(EE{i});
end


bestU=[];
bestV=[];

t = 5;
% dd1 = cache_results(@diffDist, {V,sort(diag(E)),t});
dd1 = cache_results(@diffDist, {VV{1},sort(diag(EE{1})),t});
dd2 = cache_results(@diffDist, {VV{2},sort(diag(EE{2})),t});

%%
tic
for p = 1 : length(opt.percentages)
    fprintf('%03d%% correspondences, K=%03d\n',opt.percentages(p),K);
    
    % resize D (eigenvals lists) accordingly
    Dk = D( 1:K, : );
    
    % Trim Phi and Psi
    Phi = VV{1}( :, 1:K );
    Psi = VV{2}( :, 1:K );
    
    % note that a and acc should contain the same results, this is just
    % provided to verify that our different methods for calculating
    % accuracy coincide
    acc = zeros(opt.numTests, 2);
    NMI = zeros(opt.numTests, 2);
        
    fvals = {};
    outputs = {};
    for test = 1 : opt.numTests
        fprintf('%03d%% correspondences, K=%03d, Test #%d\n',opt.percentages(p),K,test);
        if ( opt.percentages(p)~=100 || test==1 )
            if LOADCORRESP
                L = opt.previousTests.results{p}.L{test};
            else
                [L,~]=genCorrespondances_TS(opt.percentages(p),classes,dd1,dd2,VV,m,te,s,l);
            end
            %%
            
            LL = length(L);
            
            P = zeros(LL,N);
            ind = sub2ind( size(P), (1:LL)', L(:,1) );
            P(ind) = opt.alpha;
            
            Q = zeros(LL,M);
            ind = sub2ind( size(Q), (1:LL)', L(:,2));
            Q(ind) = opt.alpha;

            
            % initialization (note that here Phi=VV{1} and Psi=VV{2})
            if (initJADE)
                % this is initialization with JADE - you can also do simple
                % init with initialize (see below)
                AA = Phi\V(:,1:Z);
                BB = Psi\V(:,1:Z);

            else
                [AA, BB] = initialize( Z, P, Phi, Q, Psi);
            end
            
            X = [AA(:); BB(:)];
            u0 = Phi*AA;
            v0 = Psi*BB;
            
            optSet = optimset( ...%'PlotFcns',@optimplotfval,...
                'Display','iter',...
                'MaxIter',opt.maxIter,...
                'MaxFunEvals',1.4e6,...
                'GradObj','on',...
                'DerivativeCheck','off',...
                'GradConstr','on',...
                'Algorithm', opt.algorithm);
            
            
            [x,fval,~,output] =  fmincon(@(x)func3constr(x, P*Phi, Q*Psi, Dk, opt.lambda, opt.gamma, opt.gamma, K, Z),X,[],[],[],[],[],[],@(x)ceq_new(x,Z),optSet);
            
            % now need to derive the matrices A and B
            AA = reshape(x(1:K*Z), [K, Z]); % the potenial eigenvalues
            BB = reshape(x((K*Z+1):end), [K, Z]); % the potenial eigenvalues
            
            %finding a new basis
            u = Phi*AA;
            v = Psi*BB;
            
            % regardless of the method, we'll end up with U and V and we
            % want to calculate the accuracy of a k-means calculated in
            % those spaces
            %[~, lbl] = accuracies(u,v,clusters,classes(:)');
%% Clustering
            u = bsxfun(@times, u, 1./(sum(u.^2, 2)));
            lbl(:,1) = kmeans(u(:,1:clusters),clusters,'Replicates',100);
            
            v = bsxfun(@times, v, 1./(sum(v.^2, 2)));
            lbl(:,2) = kmeans(v(:,1:clusters),clusters,'Replicates',100);
            
            [acc(test,1), NMI(test,1), ~] = CalcMetrics(lbl(:,1), classes(:));
            [acc(test,2), NMI(test,2), ~] = CalcMetrics(lbl(:,2), classes(:));
            
            
%% store all the information you want to keep for this test
            Ls{test}=L;             % sampled correspondences
            labels{test} = lbl;     % labels found by clustering
            fvals{test} = fval;     % optimization function values
            outputs{test} = output; % optimization function outputs
            u_{test} = u;
            v_{test} = v;
            

            
        
        else
            % if 100% corresp just copy the first result along
            acc(test,:) = acc(1,:);
            NMI(test,:) = NMI(1,:);
            % probably there is no need to repeat the following
            %             Ls{test} = Ls{1};
            %             labels{test} = labels{1};
            %             fvals{test} = fvals{1};
            %             outputs{test} = outputs{1};
        end
    end

    
    
    %%
    % save all the info for this percentage
    results{p}.K = K;                  
    results{p}.L = Ls;
    results{p}.acc = acc;
    results{p}.NMI = NMI;
    results{p}.meanacc = mean(acc,1);
    results{p}.meanNMI = mean(NMI,1);
    results{p}.labels = labels;
    results{p}.fvals = fvals;
    results{p}.outputs = outputs;
    results{p}.u = u;
    results{p}.v = v;
  
    % add save code here
end
toc

if (opt.saveResults)
    save(sprintf('%s.mat',opt.saveFilename), 'opt','results','K','Z','clusters','classes');
end



function [f,g] = func3constr(x, PPhi, QKsi, D, lambda, gamma1, gamma2, K, Z)
% we'r adding terms
% to ensure the orthogonalization of the resulted bases + we'r adding terms
% to ensure the diagonalization of according Laplacian matrix

% gamma, nabla - parameters which correspond to a importance of a new basis
% to diagonalize the Laplacian matrices

% D - is a Z x 2 matrix, with each column Z eigenvalues of real Laplacians
% of according shapes

% Constrained optimization + NO demand of orthogonality of A and B (this will be implied by the constrains)


x = x(:);
N = numel(x);
KZ = K*Z;
LZ = size(PPhi,1)*Z;
ZZ = Z*Z;

% current coefficients of the basis
A = x( 1:KZ ); A = reshape( A, [K , Z] );
B = x( (KZ+1):end ); B = reshape( B, [K , Z] );

% off diagonal elements
f =     (1/LZ)*lambda*sum(sum( (PPhi*A - QKsi*B).^2 ))         + ...
    (1/ZZ)*gamma1*( sum( sum( ( A'*diag(D(1:K,1))*A - diag(D(1:Z,1)) ).^2 ) ) ) + ...
    (1/ZZ)*gamma2*( sum( sum( ( B'*diag(D(1:K,2))*B - diag(D(1:Z,2)) ).^2 ) ) );

%gradients
if nargout > 1
    gA =        (1/LZ)*lambda*(2*(PPhi)'*(PPhi)*A - 2*(PPhi)'*QKsi*B) + ...
        (1/ZZ)*4*gamma1*( diag(D(1:K,1))*A*A'*diag(D(1:K,1))*A    - diag(D(1:K,1))*A*diag(D(1:Z,1)) );
    
    gB =        (1/LZ)*lambda*(2*(QKsi)'*(QKsi)*B - 2*((QKsi)'*(PPhi*A))) + ...
        (1/ZZ)*4*gamma2*( diag(D(1:K,2))*B*B'*diag(D(1:K,2))*B    - diag(D(1:K,2))*B*diag(D(1:Z,2)) );
    
    g = [gA(:); gB(:)];
end

