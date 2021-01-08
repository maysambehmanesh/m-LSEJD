function [L,origCor] = genCorrespondances_TS(percent,classes,dd1,dd2,VV,m,tv,s,l)

% s -> t'th neighbor in tangent space
% l -> number of neighbors in tangent space



L1 = [];
L2=[];
perc = percent / 100;

% make classes a column vector
classes=classes(:);

origCor=[];     % original matching
pre_idxSize=0;

for classNum = min(classes) : max(classes)
    idx = find(classes==classNum); 
    if ~isempty(idx)
        % NOTE: rounding might eventually provide a total #points lower than
        %       the actual percentage... not a big deal but now you know it :)
        numPts = round(perc*length(idx));
        % get numPts points doing Farthest Point Sampling
        % params for fps are: starting point(s), distances, # desired points
        [randIdx,d] = fps(randi(length(idx)), dd1(idx,idx), numPts);
        
        L0=idx(randIdx(1:numPts));  % Current selected Points
        origCor=[origCor;L0];
        L1=cat(1,L1,L0);
        
        %% Using Tangent Space
        
        X1=VV{1}(:,idx);
        X2=VV{2}(:,idx);
        N=size(X1,2);
        M=size(X1,2);
        dd01=dd1(idx,idx);
        dd02=dd2(idx,idx);
        
        [neighborhood1,~]=findNeighbor(dd01,m);
        [neighborhood2,~]=findNeighbor(dd02,m);
        
        
        % Calculate Tangent Space For L
        TS1=cell(numPts,1);
        TS2=cell(numPts,1);
        
        for ii=1:numPts
            jj1=L0(ii)-pre_idxSize;
            tempX=X1(:,[jj1; neighborhood1(:,jj1)])';
            [coeff, score, latent] = pca(tempX);
            TS1{ii}=coeff(:,1:tv)';
            
            jj2=L0(ii)-pre_idxSize;
            tempX=X2(:,[jj2; neighborhood2(:,jj2)])';
            [coeff, score, latent] = pca(tempX);
            TS2{ii}=coeff(:,1:tv)';
        end
        clear tempX coeff score latent;
        
                
%         C1=ones(mm,N,numPts);
%         C2=ones(mm,M,numPts);
        B1=ones(numPts,N);
        B2=ones(numPts,M);
        for kk=1:numPts
            jj=L0(kk)-pre_idxSize;
%             C1(:,:,kk) = TS1{kk}*(X1-repmat(X1(:,jj),1,N));
            C1= TS1{kk}*(X1-repmat(X1(:,jj),1,N));
            for i=1:N
                B1(kk,i)=norm(C1(:,i));
            end
            % erase first smallest
            B1(kk,jj)=max(B1(kk,:));
            
            jj=L0(kk)-pre_idxSize;
%             C2(:,:,kk) = TS2{kk}*(X2-repmat(X2(:,jj),1,M));
            C2 = TS2{kk}*(X2-repmat(X2(:,jj),1,M));
            for i=1:M
                B2(kk,i)=norm(C2(:,i));
            end
            % erase first smallest
            B2(kk,jj)=max(B2(kk,:));  
        end
        
        
        % Find Min in Tangent Space
%         [~,minIndex1]=min(B1');
        [~,sIndex1]=sort(B1');
%         minIndex1=sIndex1(t,:);
        minIndex1=sIndex1(s:l,:); minIndex1=minIndex1(:)';
        
        minIndex1=minIndex1+pre_idxSize;  % scale index in class 2 ,3 , ....
        
%         [~,minIndex2]=min(B2');
        [~,sIndex2]=sort(B2');
%         minIndex2=sIndex2(t,:);
        minIndex2=sIndex2(s:l,:); minIndex2=minIndex2(:)';
        
        minIndex2=minIndex2+pre_idxSize;  % scale index in class 2 ,3 , ....
        
        
        
        L2=cat(1,L2,L0);
        L2=cat(1,L2,minIndex2');
        L1=cat(1,L1,minIndex1');
        
        pre_idxSize=pre_idxSize+length(idx);
    end
end

L=[L1 L2];

end