clear
close all
clc

%% Nelinearno separabilne klase

N = 500;
ro1 = rand(1,N);
teta1 = rand(1,N)*2*pi;

X = zeros(2,N);
X(1,:) = ro1.*cos(teta1);
X(2,:) = ro1.*sin(teta1);

ro2 = 1.9*rand(1,N)+2;
teta2 = rand(1,N)*2*pi;

Y(1,:) = ro2.*cos(teta2);
Y(2,:) = ro2.*sin(teta2);

figure()
hold all
scatter(X(1,:),X(2,:))
scatter(Y(1,:),Y(2,:))
title('Odbirci')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2')

%% 

X = [X Y];

pom = rand(1,2*N);
X1 = [];
X2 = [];

K = zeros(1,2*N);

for i = 1:2*N
    if pom(i) < 0.5
        K(i)=1;
    else
        K(i)=2;
    end
end

K(1:100) = 1;
K(501:600) = 2;

for i=1:2*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
    end
end
            
figure()
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
hold off
title('Slučajna klasterizacija')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2')
%%

q = ones(2*N,2)/2;

q(1:100,1) = 1; 
q(501:600,2) = 1;
q(1:100,2) = 0; 
q(501:600,1) = 0;

lmax = 1000;
l = 1;

while(l<lmax)
   Pj = zeros(1,2);
   Mj = zeros(2,1,2);
   Sj = zeros(2,2,2);
   fj = zeros(4*N,2);
   
   for j=1:2
       Pj(j)=sum(q(:,j))/(2*N);
 
       Mj(:,:,j)=zeros(2,1);
       for i=1:2*N
           Mj(:,:,j)=Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j)=Mj(:,:,j)/(2*N*Pj(j));
       
       Sj(:,:,j)=zeros(2,2);
       for i=1:2*N
           Sj(:,:,j)=Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j)=Sj(:,:,j)/(2*N*Pj(j));
       
   end
   q_pre=q;
   for i=1:2*N
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

for i=1:2*N
   [qmax,idx]=max(q_pre(i,:));
   K(i)=idx;
end

figure()
K1 = [];
K2 = [];

for i=1:2*N
    if(K(i) == 1)
        K1 = [K1, X(:,i)];
    end
    if(K(i) == 2)
        K2 = [K2, X(:,i)];
    end
    hold all
end

scatter(K1(1,:),K1(2,:))
hold all;
scatter(K2(1,:),K2(2,:))
xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])
legend('K1', 'K2')