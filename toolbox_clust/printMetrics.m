
function printMetrics(results)
    for i = 1 : length(results)
        % use this for scripts other than sparseClustering_NEG
        if any(strcmp(fieldnames(results{1}),'acc'))
            fprintf('Clustering: [%d] %s/%s (mean of %d results)\n',i,strMean(results{i}.acc),strMean(results{i}.NMI),size(results{i}.acc,1));
        else
            fprintf('[%d] %s/%s\n',i,strMean(results{i}.accuracy),strMean(results{i}.nmi));
        end
    end
    
    function [str] = strMean(data)
    % given a data matrix with more columns, it calculates mean and std of all
    % values and then prints out the best mean plusminus std
        m = mean(data,1);
        %s = std(data); 
        str = sprintf('%.1f',m(m==max(m))*100);
    end

end
        
        