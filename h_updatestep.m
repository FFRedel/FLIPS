function [h] = h_updatestep(g,h,gamma_opt,k,manual)
% Performing update step given direction and stepsize

% Remove indices that are almost zero (only applied in Linear Oracle)
if strcmp(manual,'on')
    for i=1:k
        if abs(h(i))<=0.0001 
            h(i)=0;
        end
    end
end

% h update step
h = h + gamma_opt*(g-h);

end

