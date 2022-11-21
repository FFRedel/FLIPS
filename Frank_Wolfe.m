function [H,store] = Frank_Wolfe(X,phi,epsilon,N,k,maxiter,manual,H)
% Frank-Wolfe algorithm and exact line search approach 

% Initialization
store.eta     = zeros(N,maxiter);
store.gamma   = zeros(1,maxiter);
store.sparse  = zeros(N,maxiter);
store.maxiter = maxiter;
store.h       = zeros(k,maxiter);
phi_g         = zeros(k,1);
gamma_opt     = 0;
phitrans      = phi'; 
phitphi       = phitrans*phi;

for c = 1:N
    x        = X(:,c); 
    h        = H(:,c); 
    phi_h    = phi*h;
    g        = 0*h;
    phitx    = phitrans*x;
    phitphih = phitphi*h;

    for t = 1:maxiter

        % Gradient eta computation
        [v,store,phi_h,phitphih,eta_fg,alpha] = gradient_eta(y,epsilon,c,t,store,gamma_opt,phi_h,phi_g,phitphi,phitx,phitphih,g,h);

        % Check for infeasible results 
        if alpha < 0
            disp('alpha smaller than zero ERROR, BREAK')
            break 
        end

        % Selection of descent direction
        [g,I] = g_descentdirection_FW(v); 

        % Optimal step size selection
        [gamma_opt,phi_g] = stepsize_selection(phi,g,phi_h,x,epsilon,I,alpha,eta_fg);
        store.gamma(t) = store.gamma(t)+gamma_opt; 

        % h update step
        [h] = h_updatestep(g,h,gamma_opt,k,manual);
        store.sparse(c,t) = sum(h~=0); 

        % store h
        store.h(:,t) = h;
    end

    % Store new h
    [etastar] = eta_func(x,epsilon,phi_h,0,0);
    H(:,c) = h*etastar;
end

% Total cost for all patches per iteration
if size(store.eta,1) ~= 1
    store.eta = sum(store.eta);
end

% Average stepsize for all patches per iteration
store.gamma= store.gamma/N;

end
