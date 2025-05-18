function  [M1,E,ErrM, ErrE] = ouheGP_Richardson2(alpha_u,alpha_v,tau,h)
%% 初始化定解问题

a= -10;
b= 10;
T=0.5;                %   @T 时间求解区域（暂定）
N=fix(T/tau(1));            %   @N t被分割的区间数
M=fix((b-a)/h);              %   @M x被分割的区间数
t=0:tau:T;          %   @t 时间向量
x=a:h:b;            %   @x 空间向量

D=0.2;
beta_11=1;
beta_22=1;
beta_12=0.4;
lambda=-0.3;
V=(1/2*x(1:end-1).^2).';

phi_1= 1/(sqrt(pi))*exp(-x.^2);     %  phi_1初值函数
phi_2= 1/(sqrt(pi))*exp(-x.^2);     %  phi_2 初值函数

u_temp_1=phi_1(1:M).';      %   定义u初值条件
v_temp_1=phi_2(1:M).';      %   定义v初值条件
u_temp_2=u_temp_1; v_temp_2=v_temp_1;
% 创建拉普拉斯算子矩阵A
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

% 创建算子矩阵D_1^alpha
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

% 第0层质量
M1=zeros(1,N+1);
M1(1)=h*(1/2*(abs(u_temp_1(1))^2+abs(v_temp_1(1))^2)+sum(abs(u_temp_1(2:M-1)).^2)+sum(abs(v_temp_1(2:M-1)).^2)+1/2*(abs(u_temp_1(M))^2+abs(v_temp_1(M))^2));
ErrM(1)=abs(M1(1)-M1(1));

%第0层能量
G=zeros(1,M);
G=1/2*abs(D1_u*u_temp_1).^2+1/2*abs(D1_v*v_temp_1).^2+V.*(abs(u_temp_1).^2+abs(v_temp_1).^2)...
    +D.*abs(u_temp_1).^2+1/2*beta_11*abs(u_temp_1).^4+1/2*beta_22*abs(v_temp_1).^4+beta_12*(abs(u_temp_1).^2).*(abs(v_temp_1).^2)...
    +lambda*(u_temp_1.*conj(v_temp_1)+v_temp_1.*conj(u_temp_1));
E(1)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));
ErrE(1)=abs(E(1)-E(1));

%构造第1层
u_I_1=u_temp_1-1i*tau(1)*(1/2*A_u*u_temp_1)-1i*tau(1)*(V+ones(M,1).*D+beta_11*abs(u_temp_1).^2+beta_12*abs(v_temp_1).^2).*u_temp_1...
    -1i*tau(1)*lambda*v_temp_1;
v_I_1=v_temp_1-1i*tau(1)*(1/2*A_v*v_temp_1)-1i*tau(1)*(V+beta_12*abs(u_temp_1).^2+beta_22*abs(v_temp_1).^2).*v_temp_1...
    -1i*tau(1)*lambda*u_temp_1;
uv_I_1=[u_I_1;v_I_1]; 

u_I_2=u_temp_2-1i*tau(2)*(1/2*A_u*u_temp_2)-1i*tau(2)*(V+ones(M,1).*D+beta_11*abs(u_temp_2).^2+beta_12*abs(v_temp_2).^2).*u_temp_2...
    -1i*tau(2)*lambda*v_temp_2;
v_I_2=v_temp_2-1i*tau(2)*(1/2*A_v*v_temp_2)-1i*tau(2)*(V+beta_12*abs(u_temp_2).^2+beta_22*abs(v_temp_2).^2).*v_temp_2...
    -1i*tau(2)*lambda*u_temp_2;

% 构造3层格式
I_1=ones(M,1)*(1i/(2*tau(1)));
I_2=ones(M,1)*(1i/(2*tau(2)));

 a1=I_2-1/4*V-1/4*ones(M,1)*D-1/4*beta_11*abs(u_I_2).^2-beta_12/4*(abs(v_I_2).^2);
    a2=1/2*V+1/2*ones(M,1)*D+1/2*beta_11*abs(u_I_2).^2+beta_12/2*(abs(v_I_2).^2);
    a3=I_2+1/4*V+1/4*ones(M,1)*D+1/4*beta_11*abs(u_I_2).^2+beta_12/4*(abs(v_I_2).^2);
    A1=diag(a1,0)-1/8.*A_u;
    A2=diag(a2,0)+1/4.*A_u;
    A3=diag(a3,0)+1/8.*A_u;
    u_II_2=A1\(A2*u_I_2+A3*u_temp_2+lambda*v_I_2);
    
    a4=I_2-1/4*V-1/4*beta_22*abs(v_I_2).^2-beta_12/4*(abs(u_I_2).^2);
    a5=1/2*V+1/2*beta_22*abs(v_I_2).^2+beta_12/2*(abs(u_I_2).^2);
    a6=I_2+1/4*V+1/4*beta_22*abs(v_I_2).^2+beta_12/4*(abs(u_I_2).^2);
    A4=diag(a4,0)-1/8.*A_v;
    A5=diag(a5,0)+1/4.*A_v;
    A6=diag(a6,0)+1/8.*A_v;
    v_II_2=A4\(A5*v_I_2+A6*v_temp_2+lambda*u_I_2);

        u_temp_2=u_I_2;
        u_I_2=u_II_2;   
        v_temp_2=v_I_2;
        v_I_2=v_II_2;
        uv_I_2=[u_I_2;v_I_2];   

     uv_E_I=4/3*uv_I_2-1/3*uv_I_1;
    
    %第1层之后质量
    M1(2)=h*(1/2*(abs(uv_E_I(1))^2+abs(uv_E_I(M+1))^2)+sum(abs(uv_E_I(2:M-1)).^2)+sum(abs(uv_E_I(M+2:2*M-1)).^2)+1/2*(abs(uv_E_I(M))^2+abs(uv_E_I(2*M))^2));
    ErrM(2)=abs(M1(2)-M1(1));
    %第1层之后能量
    G=1/2*abs(D1_u*uv_E_I(1:M)).^2+1/2*abs(D1_v*uv_E_I(M+1:2*M)).^2 ...
        +V.*(abs(uv_E_I(1:M)).^2+abs(uv_E_I(M+1:2*M)).^2)...
    +D.*abs(uv_E_I(1:M)).^2+1/2*beta_11*abs(uv_E_I(1:M)).^4+1/2*beta_22*abs(uv_E_I(M+1:2*M)).^4+beta_12*(abs(uv_E_I(1:M)).^2).*(abs(uv_E_I(M+1:2*M)).^2)...
    +lambda*(uv_E_I(1:M).*conj(uv_E_I(M+1:2*M))+uv_E_I(M+1:2*M).*conj(uv_E_I(1:M)));
    E(2)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));
    ErrE(2)=abs(E(2)-E(1));
    
    
for n=2:N
     a1=I_1-1/4*V-1/4*ones(M,1)*D-1/4*beta_11*abs(u_I_1).^2-beta_12/4*(abs(v_I_1).^2);
    a2=1/2*V+1/2*ones(M,1)*D+1/2*beta_11*abs(u_I_1).^2+beta_12/2*(abs(v_I_1).^2);
    a3=I_1+1/4*V+1/4*ones(M,1)*D+1/4*beta_11*abs(u_I_1).^2+beta_12/4*(abs(v_I_1).^2);
    A1=diag(a1,0)-1/8.*A_u;
    A2=diag(a2,0)+1/4.*A_u;
    A3=diag(a3,0)+1/8.*A_u;
    u_II_1=A1\(A2*u_I_1+A3*u_temp_1+lambda*v_I_1);
    
   a4=I_1-1/4*V-1/4*beta_22*abs(v_I_1).^2-beta_12/4*(abs(u_I_1).^2);
    a5=1/2*V+1/2*beta_22*abs(v_I_1).^2+beta_12/2*(abs(u_I_1).^2);
    a6=I_1+1/4*V+1/4*beta_22*abs(v_I_1).^2+beta_12/4*(abs(u_I_1).^2);
    A4=diag(a4,0)-1/8.*A_v;
    A5=diag(a5,0)+1/4.*A_v;
    A6=diag(a6,0)+1/8.*A_v;
    v_II_1=A4\(A5*v_I_1+A6*v_temp_1+lambda*u_I_1);
  
    u_temp_1=u_I_1;
    u_I_1=u_II_1;
    
    v_temp_1=v_I_1;
    v_I_1=v_II_1;
    uv_I_1=[u_I_1;v_I_1]; 
    
   for m=1:2
    a1=I_2-1/4*V-1/4*ones(M,1)*D-1/4*beta_11*abs(u_I_2).^2-beta_12/4*(abs(v_I_2).^2);
    a2=1/2*V+1/2*ones(M,1)*D+1/2*beta_11*abs(u_I_2).^2+beta_12/2*(abs(v_I_2).^2);
    a3=I_2+1/4*V+1/4*ones(M,1)*D+1/4*beta_11*abs(u_I_2).^2+beta_12/4*(abs(v_I_2).^2);
    A1=diag(a1,0)-1/8.*A_u;
    A2=diag(a2,0)+1/4.*A_u;
    A3=diag(a3,0)+1/8.*A_u;
    u_II_2=A1\(A2*u_I_2+A3*u_temp_2+lambda*v_I_2);
    
    a4=I_2-1/4*V-1/4*beta_22*abs(v_I_2).^2-beta_12/4*(abs(u_I_2).^2);
    a5=1/2*V+1/2*beta_22*abs(v_I_2).^2+beta_12/2*(abs(u_I_2).^2);
    a6=I_2+1/4*V+1/4*beta_22*abs(v_I_2).^2+beta_12/4*(abs(u_I_2).^2);
    A4=diag(a4,0)-1/8.*A_v;
    A5=diag(a5,0)+1/4.*A_v;
    A6=diag(a6,0)+1/8.*A_v;
    v_II_2=A4\(A5*v_I_2+A6*v_temp_2+lambda*u_I_2);

        u_temp_2=u_I_2;
        u_I_2=u_II_2;   
        v_temp_2=v_I_2;
        v_I_2=v_II_2;
   end 
     uv_I_2=[u_I_2;v_I_2];   
     uv_E_I=4/3*uv_I_2-1/3*uv_I_1;
    %第1层之后质量
    M1(n+1)=h*(1/2*(abs(uv_E_I(1))^2+abs(uv_E_I(M+1))^2)+sum(abs(uv_E_I(2:M-1)).^2)+sum(abs(uv_E_I(M+2:2*M-1)).^2)+1/2*(abs(uv_E_I(M))^2+abs(uv_E_I(2*M))^2));
    ErrM(n+1)=abs(M1(n+1)-M1(1));
    %第1层之后能量
    G=1/2*abs(D1_u*uv_E_I(1:M)).^2+1/2*abs(D1_v*uv_E_I(M+1:2*M)).^2 ...
        +V.*(abs(uv_E_I(1:M)).^2+abs(uv_E_I(M+1:2*M)).^2)...
    +D.*abs(uv_E_I(1:M)).^2+1/2*beta_11*abs(uv_E_I(1:M)).^4+1/2*beta_22*abs(uv_E_I(M+1:2*M)).^4+beta_12*(abs(uv_E_I(1:M)).^2).*(abs(uv_E_I(M+1:2*M)).^2)...
    +lambda*(uv_E_I(1:M).*conj(uv_E_I(M+1:2*M))+uv_E_I(M+1:2*M).*conj(uv_E_I(1:M)));
    E(n+1)=h*(1/2*G(1)+sum(G(2:M-1))+1/2*G(M));
    ErrE(n+1)=abs(E(n+1)-E(1));
  
end

    



