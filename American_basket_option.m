function American_basket_option
format long
%Input parameters and change grad, x1_points, x2_points, z_point related to different rho
n=20;       %spatial node
r1=0.03;    %r
rho=0.3;    %\rho
alpha1=0.6;
alpha2=0.4;
K=1;
[UU,t,p,Q]=multi_asset_ADMM(n,r1,rho);

[p1,s1,t1,A1]=uniform_bianhuan(15,-1,1,Q);

% The number of mesh points
n1=size(p1,2);
U=zeros(1,n1);
grad=[783,768,766,762,761,769,758,742,741,733,704,732,730,729,746,725,726,709,707,695,697,660,662,694,692,689,713,714,685,686,670,668,665,649,651,610,600,611,612,614,646,643,671,673,675,677,638,640,620,618,615,595,596,552,549,540,537,553,555,557,591,589,621,623,626,634,629,581,583,586,562,560,529,532,534,486,484,471,407,469,467,465,489,492,528,527,525,566,568,580,578,573,571,519,521,524,497,495,459,462,464,412,410,393,395,328,330,391,390,388,415,418,458,456,454,500,502,518,516,513,512,508,505,449,451,452,423,421,382,384,387,336,334,318,321,257,258,261,263,315,313,339,341,343,379,376,427,429,430,445,443,434,432,371,373,349,347,344,308,310,269,266,264,252,254,199,201,249,247,270,272,273,305,302,350,352,354,368,356,297,299,279,277,274,241,243,207,204,202,194,147,191,188,186,210,212,238,236,234,282,285,230,233,217,216,180,183,185,153,151,142,140,138,136,156,159,179,177,176,221,175,164,162,130,133,135,107,106,95,94,110,113,115,127,125,118,116,88,90,71,70,59,74,75,77,85,78,54,55,43,30,46,47,27,29,24];%when rho=0.3
%grad=[784,767,763,738,762,760,758,741,733,735,703,704,732,729,745,726,728,708,695,696,660,653,651,663,692,691,711,724,715,687,690,667,647,650,610,609,550,598,597,614,646,644,670,686,684,676,673,641,642,617,593,595,553,537,539,483,484,536,533,556,592,589,621,622,637,677,634,626,624,586,588,560,531,532,488,467,468,408,394,392,411,465,463,491,528,527,563,566,583,627,630,631,579,568,567,523,496,495,459,462,414,388,391,331,320,259,319,318,335,387,386,418,458,456,498,500,520,570,572,576,575,573,517,503,501,453,423,421,382,385,338,314,316,262,254,256,197,253,251,265,313,311,341,343,379,425,426,450,505,506,514,445,446,430,427,375,347,344,308,309,268,248,249,200,196,195,192,204,247,244,272,273,304,348,350,372,431,433,369,353,352,301,277,276,241,243,207,189,191,147,142,150,187,186,210,213,237,278,281,297,296,284,282,234,216,214,182,155,154,138,141,105,107,135,157,159,179,218,219,231,174,175,163,160,132,111,108,96,69,92,112,113,128,164,166,125,117,116,89,73,71,58,74,77,85,118,80,78,54,43,45,47,51,48,28,25];%when rho=0.5

N=size(grad,2);
for i=1:N
    k=grad(i);
    % Coefficient of basic function
    a=[p(1,t(2,k))*p(2,t(3,k))-p(1,t(3,k))*p(2,t(2,k)) p(1,t(3,k))*p(2,t(1,k))-p(1,t(1,k))*p(2,t(3,k)) p(1,t(1,k))*p(2,t(2,k))-p(1,t(2,k))*p(2,t(1,k))];
    b=[p(2,t(2,k))-p(2,t(3,k)) p(2,t(3,k))-p(2,t(1,k)) p(2,t(1,k))-p(2,t(2,k))];
    c=[p(1,t(3,k))-p(1,t(2,k)) p(1,t(1,k))-p(1,t(3,k)) p(1,t(2,k))-p(1,t(1,k))];
    s=S(p,t(:,k));
    U(1,i)=((a(1)+b(1)*p1(1,i)+c(1)*p1(2,i))*UU(t(1,k))+(a(2)+b(2)*p1(1,i)+c(2)*p1(2,i))*UU(t(2,k))+(a(3)+b(3)*p1(1,i)+c(3)*p1(2,i))*UU(t(3,k)))/(2*s);    
end

p3(1,:)=exp(Q(1,1)*p1(1,:)+Q(1,2)*p1(2,:)); 
p3(2,:)=exp(Q(2,1)*p1(1,:)+Q(2,2)*p1(2,:)); 

%%plot three-dimensional diagram
x_grid = linspace(exp(-1), exp(1), 16);   %variable x_1
y_grid = linspace(exp(-1), exp(1), 16);   %variable x_2
[x_1, x_2] = meshgrid(x_grid, y_grid); 
x = p3(1, :);  
y = p3(2, :);  
F = scatteredInterpolant(x', y', U(1,:)', 'natural'); 
U_matrix = F(x_1, x_2);                         
pp=surf(x_grid,y_grid,U_matrix,'FaceAlpha',1);
pp.EdgeColor='none';
hold on

g=max(K-alpha1*x_1-alpha2*x_2,0);
[row, col] = find(U_matrix > g); %find the corresponding node to confirm x1_points and x2_points location
%example rho=0.3
x1_points = [x_grid(6),x_grid(4),x_grid(3),x_grid(2),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1)];
x2_points = [y_grid(1),y_grid(2),y_grid(3),y_grid(4),y_grid(5),y_grid(6),y_grid(7),y_grid(8),y_grid(9),y_grid(10),y_grid(11),y_grid(12),y_grid(13),y_grid(14),y_grid(15),y_grid(16)];
z_points = [U_matrix(1,6),U_matrix(2,4),U_matrix(3, 3),U_matrix(4,2),U_matrix(5,1),U_matrix(6,1),U_matrix(7,1),U_matrix(8,1),U_matrix(9,1),U_matrix(10,1),U_matrix(11,1),U_matrix(12,1),U_matrix(13,1),U_matrix(14,1),U_matrix(15,1),U_matrix(16,1)];
%example rho=0.5
% x1_points = [x_grid(4),x_grid(3),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1),x_grid(1)];
% x2_points = [y_grid(1),y_grid(2),y_grid(3),y_grid(4),y_grid(5),y_grid(6),y_grid(7),y_grid(8),y_grid(9),y_grid(10),y_grid(11),y_grid(12),y_grid(13),y_grid(14),y_grid(15),y_grid(16)];
% z_points = [U_matrix(1,4),U_matrix(2,3),U_matrix(3, 1),U_matrix(4,1),U_matrix(5,1),U_matrix(6,1),U_matrix(7,1),U_matrix(8,1),U_matrix(9,1),U_matrix(10,1),U_matrix(11,1),U_matrix(12,1),U_matrix(13,1),U_matrix(14,1),U_matrix(15,1),U_matrix(16,1)];
plot3(x1_points,x2_points,z_points,'k-','LineWidth',3)
xlabel('$S_1$','Interpreter','latex','fontsize',18); 
ylabel('$S_2$','Interpreter','latex','fontsize',18); 
zlabel('$P$','Interpreter','latex','fontsize',18,'rotation',1) 
end


function [UU,t_new,p_new,Q]=multi_asset_ADMM(n,r1,rho)
format long;
global K alpha1 alpha2 left right xi1 xi2 eta1 eta2;
%%%%%%% common coefficient %%%%%%%%%%%%
K=1;
T=1;
q1=0.03;
q2=0.03;
sigma1=0.2;
sigma2=0.4;
alpha1=0.6;
alpha2=0.4;
gamma1=0.5*sigma1^2;
gamma2=0.5*sigma2^2;
gamma3=rho*sigma1*sigma2;
mu1=r1-q1-gamma1;
mu2=r1-q2-gamma2;
H=[gamma1,0.5*gamma3;0.5*gamma3,gamma2];
Q=orth(H); 
xi1=Q(1,1);
xi2=Q(2,1);
eta1=Q(1,2);
eta2=Q(2,2);

Value=Q'*H*Q;
nu1=Value(1,1);
nu2=Value(2,2);

lamdel1=xi1*mu1+xi2*mu2;   
lamdel2=eta1*mu1+eta2*mu2; 
aa=lamdel1/nu1;
bb=lamdel2/nu2;
wg=[1/20 1/20 1/20 9/20 2/15 2/15 2/15]';
%%%%%%%%%%%%%%%%%%%%%% matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=1.25;               %domain length
h=2*L/n;
left=-L;
right=L;
[p1,e,t,C]=uniform(n,left,right);
n1=size(p1,2);        %The number of mesh points
p=p1;
n2=0;                 %The number of mesh points on the boundary
for i=1:n1
    n2=n2+p(3,i);
end
N=size(t,2);          %The number of small triangles

A1=sparse(n1,n1);
B1=sparse(n1,n1);
C1=sparse(n1,n1);

F=zeros(n1,1);         %The right-hand of all mesh points
F_1=zeros(n1-n2,1);    %The right-hand of inner mesh points
r=p(3,:);              %The flag of inner mesh points and boundary points

M=ceil(10*T/(h));      %partitions in time
delta_tau=T/M;
tau=[0:delta_tau:T];
U=zeros(M+1,n1);       %The numerical solutions of all mesh points
U_1=zeros(M+1,n1-n2);  %The numeical solutions of inner mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assignment for A and B and Computing stiffness matrix
for i=1:N
% Numeical integral node in small triangle i
    q=[p(:,t(1,i)) p(:,t(2,i)) p(:,t(3,i)) (p(:,t(1,i))+p(:,t(2,i))+p(:,t(3,i)))/3 (p(:,t(1,i))+p(:,t(2,i)))/2 (p(:,t(1,i))+p(:,t(3,i)))/2 (p(:,t(2,i))+p(:,t(3,i)))/2];
    s=S(p,t(:,i));
    q1=aa*q(1,:)+bb*q(2,:); 
    q2=exp(q1);                %kappa every triangle seven point
% Coefficient of basic function
    a=[p(1,t(2,i))*p(2,t(3,i))-p(1,t(3,i))*p(2,t(2,i)) p(1,t(3,i))*p(2,t(1,i))-p(1,t(1,i))*p(2,t(3,i)) p(1,t(1,i))*p(2,t(2,i))-p(1,t(2,i))*p(2,t(1,i))];
    b=[p(2,t(2,i))-p(2,t(3,i)) p(2,t(3,i))-p(2,t(1,i)) p(2,t(1,i))-p(2,t(2,i))];
    c=[p(1,t(3,i))-p(1,t(2,i)) p(1,t(1,i))-p(1,t(3,i)) p(1,t(2,i))-p(1,t(1,i))];
    x=ones(size(wg'))*wg*s;
% Loop for mesh points of the i's small triangle    
    for j=1:3 
% Determining whether the mesh point is inner point or boundary point
        if r(t(j,i))==0
            for l=1:3 
% Determining whether the mesh point is inner point or boundary point
                if  r(t(l,i))==0 
                    A1(t(j,i),t(l,i))=s*q2.*phi(a(l),b(l),c(l),q,s).*phi(a(j),b(j),c(j),q,s)*wg+ A1(t(j,i),t(l,i));
                    B1(t(j,i),t(l,i))=s*q2.*phi_x(b(l),q,s).*phi_x(b(j),q,s)*wg+ B1(t(j,i),t(l,i));
                    C1(t(j,i),t(l,i))=s*q2.*phi_y(c(l),q,s).*phi_y(c(j),q,s)*wg+ C1(t(j,i),t(l,i));
                end
            end
        end
    end
end
M_globel=(1+r1*delta_tau)*A1+(nu1*B1+nu2*C1)*delta_tau;
N_globel=A1;

% Extracting the inner matrix element
M_globel(find(sum(abs(M_globel),2)==0),:)=[];
M_globel(:,find(sum(abs(M_globel),1)==0))=[]; 
M_1globel=M_globel;
N_globel(find(sum(abs(N_globel),2)==0),:)=[];
N_globel(:,find(sum(abs(N_globel),1)==0))=[]; 
N_1globel=N_globel; 
A=M_1globel;  
B1=N_1globel;
% The initial value
j=1;
 for i=1:n1
     U(1,i)=max(K-alpha1*exp(xi1*p(1,i)+eta1*p(2,i))-alpha2*exp(xi2*p(1,i)+eta2*p(2,i)),0);
     if r(i)==0
        U_1(1,j)=U(1,i);
        MM(j)=U(1,i);
        j=j+1;           
    end
 end
G=MM';
% The boundary condition
for i=1:n1
    if r(i)>0
        U(:,i)=max(K-alpha1*exp(xi1*p(1,i)+eta1*p(2,i))-alpha2*exp(xi2*p(1,i)+eta2*p(2,i)),0);
     end
end

for k=1:M 
    F=zeros(n1,1);
    F_1=zeros(n1-n2,1);
    for i=1:N  % Loop for small triangles
    % Numeical integral node in small triangle i
        q=[p(:,t(1,i)) p(:,t(2,i)) p(:,t(3,i)) (p(:,t(1,i))+p(:,t(2,i))+p(:,t(3,i)))/3 (p(:,t(1,i))+p(:,t(2,i)))/2 (p(:,t(1,i))+p(:,t(3,i)))/2 (p(:,t(2,i))+p(:,t(3,i)))/2];
        s=S(p,t(:,i));
        q1=aa*q(1,:)+bb*q(2,:);
        q2=exp(q1);
    % Coefficient of basic function
        a=[p(1,t(2,i))*p(2,t(3,i))-p(1,t(3,i))*p(2,t(2,i)) p(1,t(3,i))*p(2,t(1,i))-p(1,t(1,i))*p(2,t(3,i)) p(1,t(1,i))*p(2,t(2,i))-p(1,t(2,i))*p(2,t(1,i))];
        b=[p(2,t(2,i))-p(2,t(3,i)) p(2,t(3,i))-p(2,t(1,i)) p(2,t(1,i))-p(2,t(2,i))];
        c=[p(1,t(3,i))-p(1,t(2,i)) p(1,t(1,i))-p(1,t(3,i)) p(1,t(2,i))-p(1,t(1,i))];
    % Loop for mesh points of the i's small triangle    
        for j=1:3 
    % Determining whether the mesh point is inner point or boundary point
            if r(t(j,i))==0
                for l=1:3 
    % Determining whether the mesh point is inner point or boundary point
                    if  r(t(l,i))>0 
                    F(t(j,i))=F(t(j,i))+(U(k+1,t(l,i))-U(k,t(l,i)))*s*(q2.*phi(a(l),b(l),c(l),q,s).*phi(a(j),b(j),c(j),q,s)*wg)...
                    +delta_tau*U(k+1,t(l,i))*(nu1*q2.*phi_x(b(j),q,s).*phi_x(b(l),q,s)...
                    +nu2*q2.*phi_y(c(l),q,s).*phi_y(c(j),q,s)...
                    +r1*q2.*phi(a(l),b(l),c(l),q,s).*phi(a(j),b(j),c(j),q,s))*wg*s;
                    end
                end
            end
        end
    end              
 m=1;
    for i=1:n1
       if r(i)==0  
         F_1(m)=F(i);
         m=m+1;           
       end
    end 
end
%%%%%%%%%%%%%%%%%  main procedure 
[rows] = size(A);            %the matrix A is Coefficient matrix of U^{m}
mu=1e-4;                     %penalty parameter
P=inv(A+mu*eye(rows));       %(A+mu*I)^{-1}
err=0.1*h^(5/2);             %inexact stopping criterion

for i=1:M 
    B=-B1*U_1(i,:)'+F_1;                  %B^{m}(U^{m-1})        
    Z_cur=max(G,U_1(i,:)');
    U_1(i+1,:)=ADMM(P,B,G,Z_cur,mu,err);  %ADMM iteration 
end

j=1;
 for i=1:n1
     if r(i)==0 
         U(:,i)=U_1(:,j);
         j=j+1;
     end
 end
UU=U(M+1,:);

le=size(t,2);
t_new = zeros(4,le);
t_new(1:3,:)=t(:,:);
for i=1:le
    t_new(4,i)=1;
end
p_new = p(1:2,:);
end

function  Phi_Yd=ADMM(P,B,G1,Z_cur,mu,err) % ADMM Algorithm
    La_pre = 1*ones(size(B));  
    La_cur = zeros(size(B));   
    Z_pre = zeros(size(B));  
    k = 1;
    while  max(norm(Z_cur-Z_pre,inf),norm(La_pre-La_cur,inf))>err 
          Z_pre = Z_cur;
          La_pre = La_cur;           
          Y_pre = P*(-B-La_pre+mu*Z_pre);     %Y subproblem               
          Z_cur = max(Y_pre+1/mu*La_pre,G1);  %Z subproblem  
          La_cur = La_pre-mu*(Z_cur-Y_pre);   %Lamda subproblem  
          k = k+1;                  
    end
     Phi_Yd = Z_cur;
end

% Basis function in small triangle
function f=phi(a,b,c,q,s)
f=(0.5/s)*(a+b*q(1,:)+c*q(2,:));
end

% The derivative of basis function in small triangle
function f=phi_x(b,q,s)
f=0.5*b/s*ones(size(q(1,:)));
end
  
function f=phi_y(c,q,s)
f=0.5*c/s*ones(size(q(2,:)));
end

% Triangle area
function f=S(p,t)
M=[1 p(1,t(1)) p(2,t(1));1 p(1,t(2)) p(2,t(2));1 p(1,t(3)) p(2,t(3))];
f=1/2*det(M);
end
