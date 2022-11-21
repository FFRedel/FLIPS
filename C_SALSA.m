function [H,store] = C_SALSA(X,phi,epsilon,N,k,H,maxiter,mu1)
% C-SALSA

% Initialization
store.eta     = zeros(N,maxiter);
store.gamma   = zeros(1,maxiter);
store.sparse  = zeros(N,maxiter);
store.h       = zeros(k,maxiter);
store.eta     = zeros(N,maxiter);
store.gamma   = zeros(1,maxiter);
store.maxiter = maxiter;
phiphi        = phi'*phi;
invphiphi     = inv(eye(size(phi)) + phiphi);

tau= 1/mu1;
for c = 1:N
    
    % Initialization for every iteration
    x  = X(:,c); 
    h  = H(:,c); 
    v1 = 0*h;
    v2 = 0*x;
    d1 = 0*h;
    d2 = 0*x;
    
    psi_function = @(h,tau) soft(h,tau);
    phi_function = @(h) sum(abs(h(:))); 
    
    for t=1:maxiter
        if t ==1 
            if norm(x,2)^2<epsilon^2
                u = 0*h;
                break
            end
        end

        r = v1 + d1 + phi'*(v2 + d2);
        u = invphiphi*r;
        v1= psi_function(u-d1,tau);
        Du = phi*u;
        s = Du-d2;
        if norm(s-x,2)> epsilon
            v2= x+ epsilon*(s-x)/norm(s-x,2);
        else
            v2= x+ s-x;
        end
        d1 = d1 - u  + v1;  % Lagrange variable 1
        d2 = d2 - Du + v2;  % Lagrange variable 2
        store.eta(c,t) = phi_function(u);

    end
    H(:,c) = u;
 end

% Total cost for all patches per iteration
if size(store.eta,1) ~= 1
    store.eta = sum(store.eta);
end

% Reduce memory storage
H = sparse(H);
end
 