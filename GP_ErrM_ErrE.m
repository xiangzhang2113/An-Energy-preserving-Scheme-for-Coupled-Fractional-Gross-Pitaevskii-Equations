function  [ErrM, ErrE] = GP_ErrM_ErrE(alpha_u,alpha_v,tau,h)
% Initialize problem domain
a= -10; b= 10;       % Spatial domain [a, b]     
T=0.5;               % Final simulation time
N=fix(T/tau);        % Number of time steps
M=fix((b-a)/h);      % Number of spatial steps
t=0:tau:T;           % Time vector
x=a:h:b;             % Spatial grid

% System parameters
D=0.2;
beta_11=1;
beta_22=1;
beta_12=0.4;
lambda=-0.3;
V=(1/2*x(1:M).^2).';

% Initial conditions (Gaussian functions)
phi_1= 1/(sqrt(pi))*exp(-x.^2); 
phi_2= 1/(sqrt(pi))*exp(-x.^2);  
u_temp=phi_1(1:M).';      % Initial condition for u
v_temp=phi_2(1:M).';      % Initial condition for v
uv_temp=[u_temp;v_temp];  % Combined state vector

% Construct A (Fractional Laplacian)
r=ones(1,M);
for j=1:M
     r(j)=1/(2*M)*(abs(-M/2*(2*pi/(b-a))))^alpha_u*exp(1i*(-M/2*(2*pi/(b-a)))*((j-1)*h))+...
          1/(2*M)*(abs(M/2*(2*pi/(b-a))))^alpha_u*exp(1i*(M/2*(2*pi/(b-a)))*((j-1)*h));
    for k=-M/2+1:M/2-1
    r(j)=r(j)+1/M*(abs(k*(2*pi/(b-a))))^alpha_u*exp(1i*(k*(2*pi/(b-a)))*(j-1)*h);
    end
end
A_u=toeplitz(r);

for j=1:M
     r(j)=1/(2*M)*(abs(-M/2*(2*pi/(b-a))))^alpha_v*exp(1i*(-M/2*(2*pi/(b-a)))*((j-1)*h))+...
          1/(2*M)*(abs(M/2*(2*pi/(b-a))))^alpha_v*exp(1i*(M/2*(2*pi/(b-a)))*((j-1)*h));
    for k=-M/2+1:M/2-1
    r(j)=r(j)+1/M*(abs(k*(2*pi/(b-a))))^alpha_v*exp(1i*(k*(2*pi/(b-a)))*(j-1)*h);
    end
end
A_v=toeplitz(r);

% Construct D1  (Fractional derivative operator)
d=ones(1,M);
for j=1:M
     d(j)=1/(2*M)*(abs(-M/2*(2*pi/(b-a))))^(alpha_u/2)*exp(1i*(-M/2*(2*pi/(b-a)))*((j-1)*h))+...
          1/(2*M)*(abs(M/2*(2*pi/(b-a))))^(alpha_u/2)*exp(1i*(M/2*(2*pi/(b-a)))*((j-1)*h));
    for k=-M/2+1:M/2-1
    d(j)=d(j)+1/M*(abs(k*(2*pi/(b-a))))^(alpha_u/2)*exp(1i*(k*(2*pi/(b-a)))*(j-1)*h);
    end
end
D1_u=toeplitz(d);

for j=1:M
     d(j)=1/(2*M)*(abs(-M/2*(2*pi/(b-a))))^(alpha_v/2)*exp(1i*(-M/2*(2*pi/(b-a)))*((j-1)*h))+...
          1/(2*M)*(abs(M/2*(2*pi/(b-a))))^(alpha_v/2)*exp(1i*(M/2*(2*pi/(b-a)))*((j-1)*h));
    for k=-M/2+1:M/2-1
    d(j)=d(j)+1/M*(abs(k*(2*pi/(b-a))))^(alpha_v/2)*exp(1i*(k*(2*pi/(b-a)))*(j-1)*h);
    end
end
 
D1_v=toeplitz(d);

%First Time Step Calculation
u_I=u_temp-1i*tau*(1/2*A_u*u_temp)-1i*tau*(V+ones(M,1).*D+beta_11*abs(u_temp).^2+beta_12*abs(v_temp).^2).*u_temp...
    -1i*tau*lambda*v_temp;
v_I=v_temp-1i*tau*(1/2*A_v*v_temp)-1i*tau*(V+beta_12*abs(u_temp).^2+beta_22*abs(v_temp).^2).*v_temp...
    -1i*tau*lambda*u_temp;
uv_I=[u_I;v_I]; 

% Compute mass at t = 0 and t = 1
M1=zeros(1,N+1);
M1(1)=h*(1/2*(abs(u_temp(1))^2+abs(v_temp(1))^2)+sum(abs(u_temp(2:M-1)).^2)+sum(abs(v_temp(2:M-1)).^2)+1/2*(abs(u_temp(M))^2+abs(v_temp(M))^2));
M1(2)=h*(1/2*(abs(u_I(1))^2+abs(v_I(1))^2)+sum(abs(u_I(2:M-1)).^2)+sum(abs(v_I(2:M-1)).^2)+1/2*(abs(u_I(M))^2+abs(v_I(M))^2));
% Compute mass conservation error
ErrM(1)=abs(M1(1)-M1(1));
ErrM(2)=abs(M1(2)-M1(1));

% Compute energy at t = 0 and t = 1
G=zeros(1,M);
G=1/2*abs(D1_u*u_temp).^2+1/2*abs(D1_v*v_temp).^2+V.*(abs(u_temp).^2+abs(v_temp).^2)...
    +D.*abs(u_temp).^2+1/2*beta_11*abs(u_temp).^4+1/2*beta_22*abs(v_temp).^4+beta_12*(abs(u_temp).^2).*(abs(v_temp).^2)...
    +lambda*(u_temp.*conj(v_temp)+v_temp.*conj(u_temp));
E(1)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));

% Compute energy conservation error
ErrE(1)=abs(E(1)-E(1));

G=1/2*abs(D1_u*u_I).^2+1/2*abs(D1_v*v_I).^2+V.*(abs(u_I).^2+abs(v_I).^2)...
    +D.*abs(u_I).^2+1/2*beta_11*abs(u_I).^4+1/2*beta_22*abs(v_I).^4+beta_12*(abs(u_I).^2).*(abs(v_I).^2)...
    +lambda*(u_I.*conj(v_I)+v_I.*conj(u_I));
E(2)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));
ErrE(2)=abs(E(2)-E(1));

% Construct the third time step calculation
I=ones(M,1)*(1i/(2*tau));
for n=2:N
    a1=I-1/4*V-1/4*ones(M,1)*D-1/4*beta_11*abs(u_I).^2-beta_12/4*(abs(v_I).^2);
    a2=1/2*V+1/2*ones(M,1)*D+1/2*beta_11*abs(u_I).^2+beta_12/2*(abs(v_I).^2);
    a3=I+1/4*V+1/4*ones(M,1)*D+1/4*beta_11*abs(u_I).^2+beta_12/4*(abs(v_I).^2);
    A1=diag(a1,0)-1/8.*A_u;
    A2=diag(a2,0)+1/4.*A_u;
    A3=diag(a3,0)+1/8.*A_u;
   
    A4=diag(ones(M,1)*lambda,0);
    
    b1=I-1/4*V-1/4*beta_22*abs(v_I).^2-beta_12/4*(abs(u_I).^2);
    b2=1/2*V+1/2*beta_22*abs(v_I).^2+beta_12/2*(abs(u_I).^2);
    b3=I+1/4*V+1/4*beta_22*abs(v_I).^2+beta_12/4*(abs(u_I).^2);
    B1=diag(b1,0)-1/8.*A_v;
    B2=diag(b2,0)+1/4.*A_v;
    B3=diag(b3,0)+1/8.*A_v;
    
    AB1=[A1,A4.*(-1/4);A4.*(-1/4),B1];
    AB2=[A2,A4.*(1/2);A4.*(1/2),B2];
    AB3=[A3,A4.*(1/4);A4.*(1/4),B3];
    uv_II=AB1\( AB2*uv_I+AB3*uv_temp);
    
    % Mass and energy update
    M1(n+1)=h*(1/2*(abs(uv_II(1))^2+abs(uv_II(M+1))^2)+sum(abs(uv_II(2:M-1)).^2)+sum(abs(uv_II(M+2:2*M-1)).^2)+1/2*(abs(uv_II(M))^2+abs(uv_II(2*M))^2));
    ErrM(n+1)=abs(M1(n+1)-M1(1));
    
    G=1/2*abs(D1_u*uv_II(1:M)).^2+1/2*abs(D1_v*uv_II(M+1:2*M)).^2 ...
        +V.*(abs(uv_II(1:M)).^2+abs(uv_II(M+1:2*M)).^2)...
    +D.*abs(uv_II(1:M)).^2+1/2*beta_11*abs(uv_II(1:M)).^4+1/2*beta_22*abs(uv_II(M+1:2*M)).^4+beta_12*(abs(uv_II(1:M)).^2).*(abs(uv_II(M+1:2*M)).^2)...
    +lambda*(uv_II(1:M).*conj(uv_II(M+1:2*M))+uv_II(M+1:2*M).*conj(uv_II(1:M)));
    E(n+1)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));
    ErrE(n+1)=abs(E(n+1)-E(1));

     % Update time layers
    uv_temp=uv_I;
    uv_I=uv_II;
    u_I=uv_I(1:M); v_I=uv_I(M+1:2*M);
end
ErrM=ErrM(N+1);
ErrE=ErrE(N+1);

