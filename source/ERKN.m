function x = ERKN(ButcherArray, f, T, x0)
    % Returns the iterations of an ERK method
    % ButcherArray: Struct with the ERK's Butcher array
    % f: Function handle
    %    Vector field of ODE, i.e., x_dot = f(t,x)
    % T: Vector of time points, 1 x Nt
    % x0: Initial state, Nx x 1
    % x: ERK iterations, Nx x Nt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = ButcherArray.A;  
    b = ButcherArray.b;
    c = ButcherArray.c;
    
    Nx = size(x0, 1);
    Nt = size(T, 2);
    Nstage = size(c,1)
    
    K = zeros(Nx, Nstage);
    dT = diff(T);
    
    x = zeros(Nx, Nt);
    x(:,1) = x0;
    
    % For each x(t), calculate x(t+1) = x(t) + delta_T*SUM(bi*ki)
    % This sum can be calulated by matrix multiplication with
    % delta_T*K*b = delta_T*[k1 k2 .. kn]*b
    for nt=1:Nt-1
        % Define variables for x(t)
        xt = x(:,nt);
        tk = T(nt);
        dTk = dT(nt);
        
        K(:,1) = f(tk, xt);
        for nstage=1:Nstage-1 % Calculate the Ks
            ninterval = 1:nstage;
            t_next = tk + dTk*c(nstage+1);
            x_next = xt + dTk*K(:,ninterval)*A(nstage+1, ninterval)';
            K(:,nstage+1) = f(t_next, x_next);
        end
        x(:,nt+1) = xt + dTk*K*b; % Use K to calculate x(t)
    end
end

