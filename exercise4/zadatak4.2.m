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

%% Prikaz odbiraka

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

%% Primena algoritma klasterizacije

X = [K1; K2; K3; K4]';

pom = rand(1,4*N);

K1 = [];
K2 = [];
K3 = [];
K4 = [];

K = zeros(1,4*N);

for i = 1:4*N
    if pom(i) < 0.25
        K(i)=1;
    elseif pom(i) <0.5
        K(i)=2;
    elseif pom(i)< 0.75
        K(i)=3;
    else 
        K(i)=4;
    end
end

K(1:100) = 1;
K(501:600) = 2;
K(1001:1100) = 3;
K(1501:1600) = 4;

for i=1:4*N
    switch(K(i))
        case 1 
            K1 = [K1, X(:,i)];
        case 2 
            K2 = [K2, X(:,i)];
        case 3 
            K3 = [K3, X(:,i)];
        case 4 
            K4 = [K4, X(:,i)];
    end
end        

figure(2)
scatter(K1(1,:), K1(2,:))
hold on
scatter(K2(1,:), K2(2,:)')
scatter(K3(1,:), K3(2,:))
scatter(K4(1,:), K4(2,:))
hold off
title('Slučajna klasterizacija')
xlabel('x1')
ylabel('x2')

%%

q = ones(4*N,4)/4;

q(1:100,1) = 1; 
q(501:600,2) = 1;
q(1001:1100,3) = 1;
q(1501:1600,4) = 1;
q(1:100,2) = 0; 
q(501:600,1) = 0;
q(1001:1100,1) = 0;
q(1501:1600,1) = 0;
q(1:100,3) = 0; 
q(501:600,3) = 0;
q(1001:1100,2) = 0;
q(1501:1600,2) = 0;
q(1:100,4) = 0; 
q(501:600,4) = 0;
q(1001:1100,4) = 0;
q(1501:1600,3) = 0;

lmax = 100;
l = 1;

while(l<lmax)
   Pj = zeros(1,4);
   Mj = zeros(2,1,4);
   Sj = zeros(2,2,4);
   fj = zeros(4*N,4);
   
   for j=1:4
       Pj(j) = sum(q(:,j))/(4*N);
 
       Mj(:,:,j) = zeros(2,1);
       for i = 1:4*N
           Mj(:,:,j) = Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j) = Mj(:,:,j)/(4*N*Pj(j));
       
       Sj(:,:,j) = zeros(2,2);
       for i=1:4*N
           Sj(:,:,j) = Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j) = Sj(:,:,j)/(4*N*Pj(j));
       
   end

   q_pre = q;
   for i=1:4*N
       for j=1:4
           fj(i,j) = 1/(2*pi*det(Sj(:,:,j))^(0.5))*exp(-0.5*(X(:,i)-Mj(:,:,j))'*inv(Sj(:,:,j))*(X(:,i)-Mj(:,:,j)));
       end
       
       for j=1:4
           q(i,j) = Pj(j)*fj(i,j)/(Pj(1)*fj(i,1)+Pj(2)*fj(i,2)+Pj(3)*fj(i,3)+Pj(4)*fj(i,4));    
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

for i=1:4*N
    if(K(i) == 1)
        K1 = [K1, X(:,i)];
    end
    if(K(i) == 2)
        K2 = [K2, X(:,i)];
    end
    if(K(i) == 3)
        K3 = [K3, X(:,i)];
    end
    if(K(i) == 4)
        K4 = [K4, X(:,i)];
    end
    hold all
end

scatter(K1(1,:),K1(2,:))
hold all;
scatter(K2(1,:),K2(2,:))
scatter(K3(1,:),K3(2,:))
scatter(K4(1,:),K4(2,:))

xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])
