%Table 1. Numerical errors and convergence orders in time with h = 0.01 at final time T = 0.5. 
h = 1/100;  % Spatial step size
tau = [1/40, 1/80, 1/160, 1/320, 1/640];  % Time step sizes
alpha_u = [2, 2, 2, 1.8, 1.2];            % Fractional order for component u
alpha_v = [2, 1.8, 1.5, 2, 1.8];          % Fractional order for component v
P = length(alpha_u);
K = length(tau);

% Initialize cell arrays to store convergence results for each parameter set
rate_l2_all = cell(P,1);    % L2 norm convergence rates
rate_max_all = cell(P,1);   % Max norm convergence rates

% Loop over each (alpha_u, alpha_v) parameter pair
for p = 1:P
    fprintf('Processing alpha pair #%d: alpha_u = %.2f, alpha_v = %.2f\n', p, alpha_u(p), alpha_v(p));
    
    % Initialize containers for current parameter pair
    u_II_cell = cell(K,1);  % Store numerical solutions for u
    epsilon_cell = cell(K-1,1);   % L2 norm error between successive tau
    maxepsilon_cell = cell(K-1,1); % Max norm error between successive tau

    % Compute numerical solutions for all time steps
    [u_II_cell{1}, v_II_cell{1}] = fgp_solver(alpha_u(p), alpha_v(p), tau(1), h);
    [u_II_cell{K}, v_II_cell{K}] = fgp_solver(alpha_u(p), alpha_v(p), tau(K), h);
    % Compute L2 and max norm errors between successive time steps
    for n = 2:K-1
        [u_II_cell{n}, v_II_cell{n}] = fgp_solver(alpha_u(p), alpha_v(p), tau(n), h);
        epsilon_cell{n-1} = sqrt(h * sum( abs(u_II_cell{n-1} - u_II_cell{n}).^2 ));
        maxepsilon_cell{n-1} = max(abs(u_II_cell{n-1} - u_II_cell{n}));
    end
    epsilon_cell{K-1} = sqrt(h * sum( abs(u_II_cell{K-1} - u_II_cell{K}).^2 ));
    maxepsilon_cell{K-1} = max(abs(u_II_cell{K-1} - u_II_cell{K}));

   % Compute convergence rates using log2(epsilon_n / epsilon_{n+1})
    rate_l2 = zeros(K-2,1);
    rate_max = zeros(K-2,1);
    for n = 1:K-2
        rate_l2(n) = log2(epsilon_cell{n} / epsilon_cell{n+1});
        rate_max(n) = log2(maxepsilon_cell{n} / maxepsilon_cell{n+1});
    end

   % Store the computed convergence rates
    rate_l2_all{p} = rate_l2;
    rate_max_all{p} = rate_max;

    disp(['L2 convergence rate for p = ', num2str(p), ':']);
    disp(rate_l2);
    disp(['Max error convergence rate for p = ', num2str(p), ':']);
    disp(rate_max);
end


%Table 2 Numerical errors and convergence orders in space with tau= 0.001 at final time T = 0.5.
tau = 1/10^3;                  
h = [4, 2, 1, 1/2, 1/4];       
alpha_u = [1.2, 1.4, 1.8,1.8,2]; 
alpha_v = [1.3,1.4,1.7,2,2]; 
K = length(h);                     
P = length(alpha_u);               
a = -10;
b = 10;

rate_l2_all = cell(P,1);
rate_max_all = cell(P,1);

for p = 1:P
    fprintf('Processing alpha pair #%d: alpha_u = %.2f, alpha_v = %.2f\n', p, alpha_u(p), alpha_v(p));

    u_II_cell = cell(K+1,1);
    epsilon_cell = cell(K,1);
    maxepsilon_cell = cell(K,1);

    u_II_cell{K+1} = fgp_solver(alpha_u(p), alpha_v(p), tau, h(K)/2);
    u_II_cell{1}   = fgp_solver(alpha_u(p), alpha_v(p), tau, h(1));

    for n = 2:K
        u_II_cell{n} = fgp_solver(alpha_u(p), alpha_v(p), tau, h(n)); 

        M = fix((b - a) / h(n));
        u_fine = u_II_cell{n}(1:2:M); 
        u_coarse = u_II_cell{n-1};

        epsilon_cell{n-1} = sqrt(h(n-1) * sum( abs(u_fine - u_coarse).^2 ));
        maxepsilon_cell{n-1} = max(abs(u_fine - u_coarse));
    end

    M = fix((b - a) / (h(K)/2));
    u_fine = u_II_cell{K+1}(1:2:M);
    u_coarse = u_II_cell{K};

    epsilon_cell{K} = sqrt(h(K) * sum( abs(u_fine - u_coarse).^2 ));
    maxepsilon_cell{K} = max(abs(u_fine - u_coarse));

    rate_l2 = zeros(K-1, 1);
    rate_max = zeros(K-1, 1);
    for n = 1:K-1
        rate_l2(n) = log2(epsilon_cell{n} / epsilon_cell{n+1});
        rate_max(n) = log2(maxepsilon_cell{n} / maxepsilon_cell{n+1});
    end

    rate_l2_all{p} = rate_l2;
    rate_max_all{p} = rate_max;

    disp(['L2 convergence rate for p = ', num2str(p), ':']);
    disp(rate_l2);
    disp(['Max error convergence rate for p = ', num2str(p), ':']);
    disp(rate_max);
end


%Table 3 Table 3: Numerical errors and the temporal convergence orders of the extrapolation solution UnE,j for α1 = 1.2, α2 = 1.8 with h = 0.01 at T = 0.5.

h=1/100;
tau=[1/40,1/80,1/160,1/320,1/640,1/1280];
alpha_u=1.2;alpha_v=1.8;

K=size(tau,2);
u_II_cell=cell(K,1);

[u_II_cell{1,1},v_II_cell{1,1}]=fgp_solver(alpha_u,alpha_v,tau(1),h);
[u_II_cell{K,1},v_II_cell{K,1}]=fgp_solver(alpha_u,alpha_v,tau(K),h);

for n=2:K-1
[u_II_cell{n,1},v_II_cell{n,1}]=fgp_solver(alpha_u,alpha_v,tau(n),h);
u_E_cell{n-1,1}=4/3*u_II_cell{n,1}-1/3*u_II_cell{n-1,1};
end
u_E_cell{K-1,1}=4/3*u_II_cell{K,1}-1/3*u_II_cell{K-1,1};

for n=1:K-2
epsilon_u_E_cell{n,1}=sqrt(h*sum( (abs( u_E_cell{n,1}- u_E_cell{n+1,1} )).^2) );
maxepsilon_u_E_cell{n,1}=max ( abs( u_E_cell{n,1}- u_E_cell{n+1,1} ) );
end

for n=1:K-3
rate_l2(n,1)=log2(epsilon_u_E_cell{n,1}/epsilon_u_E_cell{n+1,1}) 
rate_max(n,1)=log2(maxepsilon_u_E_cell{n,1}/maxepsilon_u_E_cell{n+1,1})
end

%Fig.1. The evolution of Mass and Energy for different alpha
tau=1/100;
h=1/2;
alpha_u=[1.5,1.8,2];
alpha_v=[2,1.2,2];
k=size(alpha_u,2);
u_cell=cell(k,1);
M1_cell=cell(k,1);
E_cell=cell(k,1);
ErrM_cell=cell(k,1);
ErrE_cell=cell(k,1);
for n=1:k
[M1_cell{n,1},E_cell{n,1},ErrM_cell{n,1},ErrE_cell{n,1}] = GP_M_E(alpha_u(n),alpha_v(n),tau,h);
End

>> figure(1); % Fig.1.(a)
t= 0:tau:0.5;
plot(t,M1_cell{1,1},'b--o'...
,t,M1_cell{2,1},'r--x'...
,t,M1_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('M');
ylim([0.4, 1]);

figure(2);  % Fig.1.(b)
plot(t,E_cell{1,1},'b--o'...
,t,E_cell{2,1},'r--x'...
,t,E_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex')
xlabel('t');ylabel('E');
ylim([0.4, 0.6]);

%Fig.2. Absolute errors of mass and energy for different alpha
tau=1/1000;
h=1/2;
alpha_u=[1.5,1.8,2];
alpha_v=[2,1.2,2];
k=size(alpha_u,2);
u_cell=cell(k,1);
M1_cell=cell(k,1);
E_cell=cell(k,1);
ErrM_cell=cell(k,1);
ErrE_cell=cell(k,1);
for n=1:k
[M1_cell{n,1},E_cell{n,1},ErrM_cell{n,1},ErrE_cell{n,1}] = GP_M_E(alpha_u(n),alpha_v(n),tau,h);
end

figure(3); % Fig.2.(a)
t= 0:tau:0.5;
loglog(t,ErrM_cell{1,1},'b--o'...
,t,ErrM_cell{2,1},'r--x'...
,t,ErrM_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('ErrM');
%ylim([0.4, 1]);

figure(4); % Fig.2.(b)
loglog(t,ErrE_cell{1,1},'b--o'...
,t,ErrE_cell{2,1},'r--x'...
,t,ErrE_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('ErrE');
%ylim([0.4, 0.6]);

%Fig.3. The evolution of the absolute errors of mass and energy for different spatial grid sizes M
h=[20/32,20/64,20/128,20/256];
tau=1/10^5;
alpha_u=1.8; alpha_v=1.2;
k=size(h,2);
ErrM_cell=cell(k,1);
ErrE_cell=cell(k,1);
for n=1:k
[ErrM_cell{n,1},ErrE_cell{n,1}] = GP_ErrM_ErrE(alpha_u,alpha_v,tau,h(n));
end
ErrM=[ErrM_cell{1,1},ErrM_cell{2,1},ErrM_cell{3,1},ErrM_cell{4,1}];
ErrE=[ErrE_cell{1,1},ErrE_cell{2,1},ErrE_cell{3,1},ErrE_cell{4,1}];
>> figure(1); %Fig.3. 
 loglog(M,ErrM,'b--o'...)
,M,ErrE,'r--x')
 h=legend('ErrM','ErrE');
set(h,'Interpreter','latex') 
 xlabel('M');ylabel('ErrM, ErrE');




%Fig4/5
>> tau=[1/100,1/200];
h=1/2;
alpha_u=[1.5,1.8,2];
alpha_v=[2,1.2,2];
k=size(alpha_u,2);
u_cell=cell(k,1);
M1_cell=cell(k,1);
E_cell=cell(k,1);
for n=1:k
[M1_cell{n,1},E_cell{n,1},ErrM_cell{n,1}, ErrE_cell{n,1}] =GP_Richardson(alpha_u(n),alpha_v(n),tau,h);
end
>> figure(1); %Fig.4. (a)
t= 0:tau:0.5;
plot(t,M1_cell{1,1},'b--o'...
,t,M1_cell{2,1},'r--x'...
,t,M1_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('M');
ylim([0.4, 1]);

figure(2); %Fig.4. (b)
plot(t,E_cell{1,1},'b--o'...
,t,E_cell{2,1},'r--x'...
,t,E_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('E');
ylim([0.4, 0.6]);

figure(3); %Fig.5. (a)
t= 0:tau:0.5;
loglog(t,ErrM_cell{1,1},'b--o'...
,t,ErrM_cell{2,1},'r--x'...
,t,ErrM_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('ErrM');
%ylim([0.4, 1]);

figure(4); %Fig.5. (b)
loglog(t,ErrE_cell{1,1},'b--o'...
,t,ErrE_cell{2,1},'r--x'...
,t,ErrE_cell{3,1},'g--d')
h=legend('$\alpha_u=1.5$,$\alpha_v=2$','$\alpha_u=1.8$,$\alpha_v=1.2$','$\alpha_u=2$,$\alpha_v=2$');
set(h,'Interpreter','latex') 
xlabel('t');ylabel('ErrE');
%ylim([0.4, 0.6]);














