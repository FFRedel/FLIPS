function [eta] = eta_func(x,epsilon,phi_h,phi_g,gamma)
% Computation of eta (in min_h eta(h) s.t. ||h||_1<=TT )

phi_hgamma   = (phi_h + gamma*(phi_g-phi_h));
eta = (norm(x,2)^2-epsilon^2)/(x'*phi_hgamma + sqrt((x'*phi_hgamma)'*(x'*phi_hgamma)- norm(phi_hgamma,2)^2*(norm(x,2)^2-epsilon^2)));

end