function [g,I] = g_descentdirection_FW(gradient_eta)
% Calculating g(h) in descent direction d(h) = g(h) - h 
% using Frank-Wolfe method for Linear sOracle

g     = zeros(size(gradient_eta,1),1); 
[~,I] = max(abs(gradient_eta));
g(I)  = -sign(gradient_eta(I))*1; 
g     = sparse(g);

end
