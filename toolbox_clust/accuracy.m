function [acc]=accuracy(labels,classes,debug)

if (nargin < 3)
  debug = 0;
end

totsize = 0;
totmax = 0;

labs = unique(labels);
for i = 1:length(labs)
    labelsize = length(find(labels==labs(i)));
    subclust = unique(classes(labels==labs(i)));
    maxsize = 0;
    for j = 1:length(subclust)
        clustsize = length(find(classes(labels==labs(i))==subclust(j)));
        if clustsize>maxsize
           maxsize = clustsize;
        end
    end 
    labelprec = maxsize / labelsize;
    if debug
        fprintf('Label %d => maxsize = %d, labelsize = %d, precision = %.4f\n', labs(i), maxsize, labelsize, labelprec);
    end
    totsize = totsize + labelsize;
    totmax = totmax + maxsize;
end
    acc = totmax/totsize;
    if debug
        fprintf('Totmax = %d, Totsize = %d, accuracy = %.4f\n', totmax, totsize, acc);
    end
end