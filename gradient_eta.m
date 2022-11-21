function [grad_eta,store,phi_h,phitphih,eta_fg,alpha] = gradient_eta(x,epsilon,c,t,store,gamma_opt,phi_h,phi_g,phitphi,phitx,phitphih,g)
% Calculating the gradient of eta in (in min_h eta(h) s.t. ||h||_1<=TT )

% Calculating net phi(h) with updated h
phi_h   = (1-gamma_opt)*phi_h + gamma_opt*phi_g;
[eta_fg] = eta_func(x,epsilon,phi_h,0,0);
store.eta(c,t) = eta_fg;

% Calculate gradient of \eta over h 
alpha    = ((x-eta_fg*phi_h)'*x- epsilon^2)/((x-eta_fg*phi_h)'*phi_h)^2; 
phitphig = phitphi*g;
grad_eta = -phitx*alpha + alpha *eta_fg*((1-gamma_opt)*phitphih + gamma_opt*phitphig); 
phitphih = (1-gamma_opt)*phitphih + gamma_opt*phitphig;

end

