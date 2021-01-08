function [AC, nmi_value, error_cnt] = CalcMetrics(label, result)

result = bestMap(label, result);
error_cnt = sum(label ~= result);

% The following is the naive accuracy implementation that was provided
% with the CalcMetrics original code. It has been substituted with a 
% more precise one now, but left here as commented code to allow one to
% replicate the results in the currently submitted PAMI paper
% AC = length(find(label == result))/length(label);

AC = accuracy(label,result);

nmi_value = nmi(label, result);

