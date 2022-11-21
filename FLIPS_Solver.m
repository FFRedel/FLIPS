function [H,store] = FLIPS_Solver(X,phi,epsilon,N,k,H,maxiter,TT,betainv,theta,manual,H0_LLL)
% Quadratic Oracle

% Initialization for storing variables
phitphi         = phi'*phi;
store.eta       = zeros(N,maxiter);
store.gamma     = zeros(1,maxiter);
store.sparse    = zeros(N,maxiter);
store.eta       = zeros(N,maxiter);
store.gamma     = zeros(1,maxiter);
store.h         = zeros(k,maxiter);
store.maxiter   = maxiter;

for c = 1:N

    % Initialization for every iteration
    x         = X(:,c); 
    h         = H(:,c);
    h_notnorm = H0_LLL(:,c);
    phi_h     = phi*h;
    phi_g     = 0*phi_h;
    gamma_opt = 0;
    dt_old    = 0;
    g         = 0*h;
    phity     = phi'*x;
    phitphih  = phitphi*h;

    for t = 1:maxiter

        if t ==1 
            if norm(x,2)^2<=epsilon^2
               h = 0*h;
               break
            elseif norm(x-phi*h_notnorm,2)>=epsilon
                epsilon= norm(x-phi*h_notnorm,2)*1.1;
            end
        end

        % Gradient eta computation
        [grad_eta,store,phi_h,phitphih,eta_fg,alpha] = gradient_eta(x,epsilon,c,t,store,gamma_opt,phi_h,phi_g,phitphi,phity,phitphih,g);

        % oracle h update step
        dt    = theta * dt_old + grad_eta;
        dt_old= dt;
        h_out = h - betainv* dt;

        % Projection step
        g = ProjectOntoL1Ball(h_out, TT);
        g = sparse(g);

        % Optimal stepsize selection
        I = [];
        [gamma_opt,phi_g] = stepsize_selection(phi,g,phi_h,x,epsilon,I,alpha,eta_fg);
        store.gamma(t) = store.gamma(t)+gamma_opt;

        % Actual h update step
        [h] = h_updatestep(g,h,gamma_opt,k,manual);
        h   = sparse(h);

        % store h
        [etacost]          = eta_func(x,epsilon,phi_h,0,0);
        store.h(:,t)       = h*etacost;

    end

    % Store new h
    H(:,c) = h*etacost;
    
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

