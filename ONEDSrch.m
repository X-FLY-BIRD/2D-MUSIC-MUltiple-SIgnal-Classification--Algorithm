function out = ONEDSrch(spec, estN)
% minimum point searching for given spectrum
% input: spec describes spectral function. First column is the varying 
% angle, second column is corresponding specrel value
% input: estN is the estimated number of sources
% output: location of minimum point

angle = spec(:, 1);
temp = zeros(3, 1);
minP = [];
for i = 1:length(angle)
    temp(1, 1) = temp(2, 1);
    temp(2, 1) = temp(3, 1);
    temp(3, 1) = spec(i, 2);
    if temp(2, 1) ~= 0 && temp(2, 1) < temp(1, 1) && temp(2, 1) < temp(3, 1)
        % find minimum point
        pair = [angle(i-1), temp(2, 1)];
        minP = [minP; pair];
    end
end 

if size(minP, 1) < estN
    out = minP;
else
    [~, idx] = sort(minP(:, 2));
    out = minP(idx(1:estN), :);
end
