clear
close all
clc

%% Generisanje odbiraka

N = 500;

M1 = [2;4];
M2 = [15;4];
M3 = [10;10];
M4 = [0;10];

S1 = [0.8 0.2; 0.2 0.8];
S2 = [1 -0.2; -0.2 1];
S3 = [0.7 -0.2; -0.2 0.7];
S4 = [0.8 0.5; 0.5 0.8];

K1 = mvnrnd(M1,S1,N);
K2 = mvnrnd(M2,S2,N);
K3 = mvnrnd(M3,S3,N);
K4 = mvnrnd(M4,S4,N);

%% Pretpostavka da postoje 2 klase

X=[K1; K2; K3; K4]';
pom = rand(1,4*N);

X1 = [];
X2 = [];

K = zeros(1,4*N);

for i = 1:4*N
    if pom(i) < 0.5
        K(i)=1;
    else
        K(i)=2;
    end
end

K(1:100) = 1;
K(1001:1100) = 2;

for i=1:4*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
    end
end
            
figure(2)
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
hold off
title('Slučajna klasterizacija')
xlabel('x1')
ylabel('x2')

%%

q = ones(4*N,2)/2;

q(1:100,1) = 1; 
q(1:100,2) = 0; 
q(1001:1100,1) = 0;
q(1001:1100,2) = 1;

lmax=100;
l=1;

while(l<lmax)
   Pj=zeros(1,2);
   Mj=zeros(2,1,2);
   Sj=zeros(2,2,2);
   fj=zeros(4*N,2);
   
   for j=1:2
       Pj(j)=sum(q(:,j))/(4*N);
 
       Mj(:,:,j)=zeros(2,1);
       for i=1:4*N
           Mj(:,:,j)=Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j)=Mj(:,:,j)/(4*N*Pj(j));
       
 
       Sj(:,:,j)=zeros(2,2);
       for i=1:4*N
           Sj(:,:,j)=Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j)=Sj(:,:,j)/(4*N*Pj(j));
       
   end
   q_pre=q;
   for i=1:4*N
       for j=1:2
           fj(i,j)=1/(2*pi*det(Sj(:,:,j))^(0.5))*exp(-0.5*(X(:,i)-Mj(:,:,j))'*inv(Sj(:,:,j))*(X(:,i)-Mj(:,:,j)));
       end
       
       for j=1:2
           q(i,j)=Pj(j)*fj(i,j)/(Pj(1)*fj(i,1)+Pj(2)*fj(i,2));    
       end
       
   end
   if max(max(q-q_pre))>10^(-3)
       l=l+1;
   else
       break;
   end
  
end

for i=1:4*N
   [qmax,idx]=max(q_pre(i,:));
   K(i)=idx;
end

K1_curr = [];
K2_curr = [];

figure(3)
for i=1:4*N
    if(K(i)==1)
        K1_curr = [K1_curr, X(:,i)];
    end
    if(K(i)==2)
        K2_curr = [K2_curr, X(:,i)];
    end
    hold all
end

scatter(K1_curr(1,:),K1_curr(2,:))
hold all;
scatter(K2_curr(1,:),K2_curr(2,:))

xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])

%%

clear
close all
clc

N = 500;

M1 = [2;4];
M2 = [15;4];
M3 = [10;10];
M4 = [0;10];

S1 = [0.8 0.2; 0.2 0.8];
S2 = [1 -0.2; -0.2 1];
S3 = [0.7 -0.2; -0.2 0.7];
S4 = [0.8 0.5; 0.5 0.8];

K1 = mvnrnd(M1,S1,N);
K2 = mvnrnd(M2,S2,N);
K3 = mvnrnd(M3,S3,N);
K4 = mvnrnd(M4,S4,N);

figure(1)
hold all
scatter(K1(:,1),K1(:,2));
scatter(K2(:,1),K2(:,2));
scatter(K3(:,1),K3(:,2));
scatter(K4(:,1),K4(:,2));

hold off
xlabel('x1');
ylabel('x2');
legend('K1','K2','K3','K4');
title('Odbirci klasa');

%% Pretpostavka da imamo 3 klase

X=[K1; K2; K3; K4]';

pom = rand(1,4*N);

X1 = [];
X2=[];
X3=[];

K=zeros(1,4*N);

for i = 1:4*N
    if pom(i) < 0.33
        K(i)=1;
    elseif pom(i) <0.566
        K(i)=2;
    else
        K(i)=3;
    end
end

K(1:100)=1;
K(701:800)=2;
K(1500:1501)=3;

for i=1:4*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
        case 3 
            X3 = [X3, X(:,i)];
    end
end
            
figure(2)
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
scatter(X3(1,:), X3(2,:))
hold off
title('Slučajna klasifikacija')
xlabel('x1')
ylabel('x2')

%%

q=ones(4*N,3)/3;

q(1:100,1)=1; 
q(701:800,2)=1;
q(1501:1600,3)=1;
q(701:800,1)=0;
q(1501:1600,1)=0;
q(1501:1600,2)=0;
q(1:100,2)=0; 
q(1:100,3)=0; 
q(701:800,3)=0;

lmax=100;
l=1;

while(l<lmax)
   Pj=zeros(1,3);
   Mj=zeros(2,1,3);
   Sj=zeros(2,2,3);
   fj=zeros(4*N,3);
   
   for j=1:3
       Pj(j)=sum(q(:,j))/(4*N);
 
       Mj(:,:,j)=zeros(2,1);
       for i=1:4*N
           Mj(:,:,j)=Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j)=Mj(:,:,j)/(4*N*Pj(j));
       
 
       Sj(:,:,j)=zeros(2,2);
       for i=1:4*N
           Sj(:,:,j)=Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j)=Sj(:,:,j)/(4*N*Pj(j));
       
   end
   q_pre=q;
   for i=1:4*N
       for j=1:3
           fj(i,j)=1/(2*pi*det(Sj(:,:,j))^(0.5))*exp(-0.5*(X(:,i)-Mj(:,:,j))'*inv(Sj(:,:,j))*(X(:,i)-Mj(:,:,j)));
       end
       
       for j=1:3
           q(i,j)=Pj(j)*fj(i,j)/(Pj(1)*fj(i,1)+Pj(2)*fj(i,2)+Pj(3)*fj(i,3));    
       end
       
   end
   if max(max(q-q_pre))>10^(-3)
       l=l+1;
   else
       break;
   end
  
end

for i=1:4*N
   [qmax,idx]=max(q_pre(i,:));
   K(i)=idx;
end

figure(3)
for i=1:4*N
    if(K(i)==1)
        scatter(X(1,i),X(2,i));
    end
    if(K(i)==2)
        scatter(X(1,i),X(2,i));
    end
    if(K(i)==3)
        scatter(X(1,i),X(2,i));
    end

    hold all
end
xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])

%%

clear
close all
clc

N = 500;

M1 = [2;4];
M2 = [15;4]; 
M3 = [10;10]; 
M4 = [0;10];

S1 = [0.8 0.2; 0.2 0.8];
S2 = [1 -0.2; -0.2 1];
S3 = [0.7 -0.2; -0.2 0.7];
S4 = [0.8 0.5; 0.5 0.8];

K1 = mvnrnd(M1,S1,N);
K2 = mvnrnd(M2,S2,N);
K3 = mvnrnd(M3,S3,N);
K4 = mvnrnd(M4,S4,N);

figure(1)
hold all
scatter(K1(:,1),K1(:,2));
scatter(K2(:,1),K2(:,2));
scatter(K3(:,1),K3(:,2));
scatter(K4(:,1),K4(:,2));

hold off
xlabel('x1');
ylabel('x2');
legend('K1','K2','K3','K4');
title('Odbirci klasa');

%% Pretpostavka da postoji 5 klasa

X=[K1; K2; K3; K4]';

pom = rand(1,4*N);

X1 = [];
X2 = [];
X3 = [];
X4 = [];
X5 = [];

K = zeros(1,4*N);

for i = 1:4*N
    if pom(i) < 0.2
        K(i)=1;
    elseif pom(i) <0.4
        K(i)=2;
    elseif pom(i)< 0.6
        K(i)=3;
    elseif pom(i)<0.8
        K(i)=4;
    else
        K(i)=5;
    end
end

K(1:100)=1;
K(401:500)=2;
K(801:900)=3;
K(1201:1200)=4;
K(1601:1700)=5;

for i=1:4*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
        case 3 
            X3 = [X3, X(:,i)];
        case 4 
            X4 = [X4, X(:,i)];
        case 5
            X5 = [X5, X(:,i)];
    end
end
            
figure()
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
scatter(X3(1,:), X3(2,:))
scatter(X4(1,:), X4(2,:))
scatter(X5(1,:),X5(2,:));
hold off
title('Slučajna klasifikacija')
xlabel('x1')
ylabel('x2')

%%

q = ones(4*N,4)/4;

q(1:100,1) = 1; 
q(401:500,2) = 1;
q(801:900,3) = 1;
q(1201:1300,4) = 1;
q(1601:1700,5) = 1;
q(401:500,1) = 0;
q(801:900,1) = 0;
q(1201:1300,1) = 0;
q(1601:1700,1) = 0;
q(1:100,2) = 0; 
q(801:900,2) = 0;
q(1201:1300,2) = 0;
q(1601:1700,2) = 0;
q(1:100,3) = 0; 
q(401:500,3) = 0;
q(1201:1300,3) = 0;
q(1601:1700,3) = 0;
q(1:100,4) = 0; 
q(401:500,4) = 0;
q(801:900,4) = 0;
q(1601:1700,4) = 0;
q(1:100,5) = 0; 
q(401:500,5) = 0;
q(801:900,5) = 0;
q(1201:1300,5) = 0;

lmax = 100;
l = 1;

while(l<lmax)
   Pj = zeros(1,5);
   Mj = zeros(2,1,5);
   Sj = zeros(2,2,5);
   fj = zeros(4*N,5);
   
   for j=1:5
       Pj(j) = sum(q(:,j))/(4*N);
 
       Mj(:,:,j) = zeros(2,1);
       for i=1:4*N
           Mj(:,:,j) = Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j) = Mj(:,:,j)/(4*N*Pj(j));
 
       Sj(:,:,j)=zeros(2,2);
       for i=1:4*N
           Sj(:,:,j) = Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j) = Sj(:,:,j)/(4*N*Pj(j));
       
   end

   q_pre = q;

   for i=1:4*N
       for j=1:5
           fj(i,j) = 1/(2*pi*det(Sj(:,:,j))^(0.5))*exp(-0.5*(X(:,i)-Mj(:,:,j))'*inv(Sj(:,:,j))*(X(:,i)-Mj(:,:,j)));
       end
       
       for j=1:5
           q(i,j) = Pj(j)*fj(i,j)/(Pj(1)*fj(i,1)+Pj(2)*fj(i,2)+Pj(3)*fj(i,3)+Pj(4)*fj(i,4)+Pj(5)*fj(i,5));    
       end
       
   end
   if max(max(q-q_pre))>10^(-3)
       l = l+1;
   else
       break;
   end
  
end

for i=1:4*N
   [qmax,idx] = max(q_pre(i,:));
   K(i) = idx;
end

figure(3)

K1 = [];
K2 = [];
K3 = [];
K4 = [];
K5 = [];

for i=1:4*N
    if(K(i) == 1)
        K1 = [K1, X(:,i)];
        %scatter(X(1,i),X(2,i));
    end
    if(K(i) == 2)
        %scatter(X(1,i),X(2,i));
        K2 = [K2, X(:,i)];
    end
    if(K(i) == 3)
        %scatter(X(1,i),X(2,i));
        K3 = [K3, X(:,i)];
    end
    if(K(i) == 4)
        %scatter(X(1,i),X(2,i));
        K4 = [K4, X(:,i)];
    end
    if(K(i) == 5)
        %scatter(X(1,i),X(2,i));
        K5 = [K5, X(:,i)];
    end
    hold all
end

scatter(K1(1,:),K1(2,:))
hold all;
scatter(K2(1,:),K2(2,:))
scatter(K3(1,:),K3(2,:))
scatter(K4(1,:),K4(2,:))
scatter(K5(1,:),K5(2,:))
legend('K1','K2','K3','K4','K5');

xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])