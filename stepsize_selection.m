function [gamma_opt,phi_g] = stepsize_selection(phi,g,phi_h,x,epsilon,I,alpha0,eta_0)
% Computing optimal stepsize with exact line search for min-max form of:
% min_h ||h||_1
% s.t.  ||phi(h) - x ||< epsilon
xe = norm(x,2)^2-epsilon^2;

% Calculation of phi_g and phi_gh
if isempty(I)
    phi_g = phi*g;
else
    phi_g = sign(g(I))*phi(:,I);
end
phi_gh = phi_g-phi_h;

% Computing gamma_max
aa = (x'*phi_gh)^2 - xe*norm(phi_gh,2)^2;
if aa>=0
    gamma_max = 1;
else 
    bb = 2*((x'*phi_h)*(x'*phi_gh) - xe*phi_h'*phi_gh );
    cc = abs(x'*phi_h)^2 - xe*norm(phi_h,2)^2;
    gamma_max(1) = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa); 
    gamma_max(2) = (-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
    gamma_max    = max(gamma_max);
    if gamma_max>1
        gamma_max=1;
    end
end 

% Selecting optimal gamma
lambda  = alpha0 * (x - eta_0*phi_h);
deriv_0 = - lambda'*phi_gh;

eta_m   = eta_func(x,epsilon,phi_h,phi_g,gamma_max);
alpha   = ((x-eta_m*phi_h)'*x- epsilon^2)/((x-eta_m*phi_h)'*phi_h)^2; 
lambda  = alpha * (x - eta_m*(phi_h+gamma_max*phi_gh));
deriv_m = - lambda'*phi_gh;

if deriv_0>0
    gamma_opt = 0; % since h is optimal
elseif deriv_m<=0 && alpha>0
    gamma_opt = gamma_max;
else
    a  = aa;
    b  = bb;
    c  = (2*x'*phi_gh*x'*phi_h*phi_h'*phi_gh - norm(phi_h,2)^2*norm(x'*phi_gh,2)^2 - xe*norm(phi_h'*phi_gh,2)^2) /norm(phi_gh,2)^2;
    gamma_alg(1) = (-b + sqrt(b^2-4*a*c))/2/a;
    eta_check(1) = eta_func(x,epsilon,phi_h,phi_g,gamma_alg(1));
    gamma_alg(2) = (-b - sqrt(b^2-4*a*c))/2/a;
    eta_check(2) = eta_func(x,epsilon,phi_h,phi_g,gamma_alg(2));
    [~,idxmineta]=min(eta_check);
    gamma_opt = gamma_alg(idxmineta);
end

end
