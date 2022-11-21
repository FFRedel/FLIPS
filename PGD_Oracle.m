function [H,store] = PGD_Oracle(X,phi,epsilon,maxiter,k,mu,manual,TT,gamma0,H,beta,N)
% Accelerated Projected Gradient Descent Oracle

% Initialization
phi_trans       = phi'; 
phitphi         = phi_trans*phi;
store.eta       = zeros(N,maxiter);
store.gamma     = zeros(1,maxiter);
store.eta       = zeros(N,maxiter);
store.gamma     = zeros(1,maxiter);
store.h         = zeros(k,maxiter);
store.maxiter   = maxiter;

for c = 1:N

    % Initialization for every iteration
    x         = X(:,c); 
    h         = H(:,c);
    phi_h     = phi*h;
    phi_g     = 0*phi_h;
    gamma_opt = 0;
    dt_old    = 0;
    g         = 0*h;
    phitx     = phi_trans*x;
    phitphih  = phitphi*h;
    eps       = epsilon(c);

    for t = 1:maxiter
        % Gradient eta computation
        [grad_eta,store,phi_h,alpha,phitphih] = gradient_eta(x,eps,c,t,store,gamma_opt,phi_h,phi_g,phitphi,phitx,phitphih,g);

        % Check for infeasible results 
         if alpha < 0
            disp(['alpha smaller than zero ERROR, BREAK, on iter ',num2str(t)])
            break 
         end

        % Gradient h update step
        dt    = beta * dt_old + grad_eta; % TO BE CHANGED
        dt_old= dt;
        h_out = h - gamma0* dt;

        % Projection step
        g = ProjectOntoL1Ball(h_out, TT);
        g = sparse(g);

        % Optimal stepsize selection
        I = [];
        [gamma_opt,phi_g] = stepsize_selection(phi,g,phi_h,x,eps,mu,I);
        store.gamma(t) = store.gamma(t)+gamma_opt;

        % Actual h update step
        [h] = h_updatestep(g,h,gamma_opt,k,manual);
        h   = sparse(h);

        % store h
        store.h(:,t) = h;
    end

    % Store new h
    [etastar] = eta_func(x,eps,phi_h,0,0);
    H(:,c) = h*etastar;
end

% Average stepsize over data samples
store.gamma= store.gamma/N;

% Reduce memory storage
H = sparse(H);

% Total cost for all patches per iteration
if size(store.eta,1) ~= 1
    store.eta = sum(store.eta);
end
 
end

