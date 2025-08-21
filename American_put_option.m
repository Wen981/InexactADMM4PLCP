function American_put_option
    format long
    %Input parameters
    T_hat = 1;    %time upper bound T^{\hat}        
    K = 1;        %strike price K
    r = 0.02;     %risk-free interest rate
    q = 0.025;    %dividend
    sigma = 0.3;  %volatility
    T = T_hat*(sigma^2)/2;    
    alpha = (r-q)/(sigma^2)-1/2;     %alpha
    beta  = alpha^2+2*r/(sigma^2);   %beta
    
    %compute domain length 
    L1 = -2.5*T*alpha+0.5*sqrt(25*T^2*(alpha^2)-20*T*log(1e-9/sqrt(5*K)));  
    X=(-r+q+(sigma^2)/2-sqrt((-r+q+(sigma^2)/2)^2+2*r*sigma^2))/(-r+q-(sigma^2)/2-sqrt((-r+q+(sigma^2)/2)^2+2*r*sigma^2));   
    if (-log(X)>L1)
        L = -log(X);
    else
        L = L1;
    end   
    N = 300;                 %space division number
    M = ceil(10*T/(2*L/N));  %tau=0.1h %time division number    
    dtao = T/M;     %temporal step size, where dt is dtao
    x = linspace(-L,L,N+1); 
    tao = linspace(0,T,M+1);
    t = T_hat-2*tao./sigma^2;%位置
    t = t(M+1:-1:1);
    dx = 2*L/N;
    s = K*exp(x);


    [P_ori,b_ref] = solve_Yd(M,N,dx,dtao,K,tao,x,r,q,L,alpha,beta,s);
    image_3d(P_ori,b_ref,s,M,N,t,K);
end


function [P_ori,b_ref] = solve_Yd(M,N,dx,dtao,K,tao,x,r,q,L,alpha,beta,s)
    %give parameters and matrix
    format long
    h=dx;
    %% Stiffness matrix
    phi = diag((2/dx)*ones(1,N-1))+diag((-1/dx)*ones(1,N-2),1)+diag((-1/dx)*ones(1,N-2),-1);
    psi = diag((dx/6)*ones(1,N-2),1)+diag((2*dx/3)*ones(1,N-1))+diag((dx/6)*ones(1,N-2),-1);
    phi = sparse(phi);
    psi = sparse(psi);
    A = dtao*phi+psi;    
   
    %initialize Y_d
    U = zeros(M+1,N+1);
    b_ref =ones(M+1,1)*min(1*r/q,1);  
   
    %Attach initial boundary values to the matrix Y_d and G0
    U(1,:) = exp(alpha*x+beta*tao(1)).*max(1-exp(x),0);
    U(:,N+1) = 0;
    U(:,1) = exp(beta*tao'+alpha*(-L)).*max(1-exp(-L),0);
    G0 = (exp(alpha*x(2:N)).*max(1-exp(x(2:N)),0))';
    mu=1e-4;
    P = inv(A+mu*eye(N-1)); %(A+\mu*I)^{-1} 
    err = (dx)^(5/2); %stopping criterion

    for i = 2:M+1
        H_m = zeros(N-1,1);
        H_m(1) = (h/6-dtao/h)*exp(alpha*(-L)+beta*tao(i))*max(K-exp(-L),0)-(h/6)*exp(alpha*(-L)+beta*tao(i-1))*max(K-exp(-L),0);%H_m(1) H_m(N-1)=0       
        B = -psi*U(i-1,2:N)'+H_m; 
        G1 = exp(beta*tao(i))*G0;                   
        Z_cur = max(G1,U(i-1,2:N)');
        [Phi_Yd,b_1] = ADMM(P,B,G1,s,Z_cur,err,mu); 
        b_ref(i) = b_1;
        U(i,2:N) = Phi_Yd(1:N-1);         
    end
    b_ref = b_ref(M+1:-1:1);  %optimal exercise boundary
    P_ori = K*exp(-alpha*ones(M+1,1)*x-beta*(tao')*ones(1,N+1)).*U;  %transform back option value P  
end


function  [Phi_Yd,b_1] = ADMM(P,B,G1,s,Z_cur,err,mu) %ADMM Algorithm
    La_pre = 1*ones(size(B)); 
    La_cur = zeros(size(B));   
    Z_pre = zeros(size(B)); 
    k = 1;
    while  max(norm( Z_cur-Z_pre,inf),norm(La_pre-La_cur,inf))>err         
          Z_pre = Z_cur;
          La_pre = La_cur;        
          Y_pre = P*(-B-La_pre+mu*Z_pre);    %Y subproblem               
          Z_cur = max(Y_pre+1/mu*La_pre,G1); %Z subproblem
          La_cur = La_pre-mu*(Z_cur-Y_pre);  %Lamda subproblem
          k = k+1;                
    end
     Phi_Yd = Z_cur;
     b_cur = find(Phi_Yd-G1>0, 1);
     b_1 = s(b_cur);    
end


function [tt,S] = MeshGeneration(t,s)  
     [tt,S] = meshgrid(t,s'); %Two dimensional mesh generation    
end


function image_3d(P_ori,b_ref,s,M,N,t,K)  %plot images
    v1 = P_ori(M+1,:); %option value in initial moment    
    figure
    plot_3D(P_ori,t,s,K,b_ref)  %three-dimensional graphs of P(S,t) values  
    hold on

    figure
    plot_v0(s,v1,K)%two-dimensional graphs of initial values P(S,0) values 
end


function plot_3D(V,t,s,K,b_1)
    P1 = flipud(V);
    [tt,S] = MeshGeneration(t,s);
    pp=surf(tt,S,P1'); 
    pp.EdgeColor='none';
    hold on
    plot3(tt,b_1,K-b_1,'k-','LineWidth',3)
    set(gca,'YDir','reverse')
    axis([0 1 0 5 0 1])
    title('$option \ value$','Interpreter','latex','fontsize',18);
    xlabel('$t$','Interpreter','latex','fontsize',18)
    ylabel('$S$','Interpreter','latex','fontsize',18)
    zlabel('$P$','Interpreter','latex','fontsize',18,'rotation',1)    
end
    
function plot_v0(s,v1,K)
    plot(s,max(K-s,0),'c--')   %payoff function
    hold on
    plot(s,v1,'bo','MarkerSize',3)
    hold on
    title('$option \ value \ at \ t=0$','Interpreter','latex','fontsize',18);
    xlabel('$S$','Interpreter','latex','fontsize',18)
    ylabel('$P$','Interpreter','latex','fontsize',18,'rotation',0)
    legend('max(K-S,0)','InADMM-P')
    axis([0 3 0 1]) %[0 3 0 1]
end
