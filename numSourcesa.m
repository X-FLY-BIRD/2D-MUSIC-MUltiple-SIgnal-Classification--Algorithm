function out = numSourcesa(D)
    % D is the diagonal eigen value matrix of the matrix of empirical
    % covariance of the antenna data
    num = NaN;
    for i = 1:size(D, 1)-2
        
        if (D(i+1, i+1)-D(i, i))/(D(i+2, i+2)-D(i+1, i+1))<=0.1 && ...
            D(i+1, i+1)/D(i+2, i+2)<=0.2

            num = size(D, 1)-i-1;
            break;
        end
        
    end
    
    out = num;
end