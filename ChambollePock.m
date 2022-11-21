function [H,store] = ChambollePock(X,phi,epsilon,N,k,H,maxiter,tau,sigma,theta)
% Chambolle-Pock algorithm 

% Initialization 
phitrans      = phi';
store.eta     = zeros(N,maxiter);
store.gamma   = zeros(1,maxiter);
store.h       = zeros(k,maxiter);
store.sparse  = zeros(N,maxiter);
store.maxiter = maxiter;

for c = 1:N

    % Initialization for every iteration
    x         = X(:,c); 
    h         = H(:,c); 
    hplus     = 0*h; 
    u         = 0*h;

    for t=1:maxiter
    
        if t == 1 
            if norm(x,2)^2<epsilon^2
                h = 0*h;
                break
            end
        end

        h_bar = h - tau*phitrans*u;
        for j=1:k
            hplus(j) = soft(h_bar(j),tau);
        end
        h = hplus + theta*(hplus - h);

        phi_h = phi*h;
        u_bar = u + sigma*phi_h - sigma*x;
        u = (1 - epsilon*sigma/norm(u_bar,2))*u_bar;

        cost = (phi_h)'*u + norm(h,1) - ((x'*u) + epsilon*norm(u,2));
        store.eta(c,t) = cost;

    end

    H(:,c) = sparse(h);
end

% Total cost for all patches per iteration
if size(store.eta,1) ~= 1
    store.eta = sum(store.eta);
end

end
