function x_soft = soft(x_orig,tau)
% Soft thresholding algorithm 

if sum(abs(tau(:)))~=0
    threshold = abs(x_orig) - tau;
    x_soft = max(threshold, 0);
    x_soft = x_soft./(x_soft + tau) .* x_orig;
else
    x_soft = x_orig;
end



