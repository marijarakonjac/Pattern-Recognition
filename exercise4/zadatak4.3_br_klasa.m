clear
close all
clc

%% Generisanje nelinearno separabilnih klasa

N = 500;
ro1 = rand(1,N);
teta1 = rand(1,N)*2*pi;

X = zeros(2,N);
X(1,:) = ro1.*cos(teta1);
X(2,:) = ro1.*sin(teta1);

ro2 = rand(1,N)+2;
teta2 = rand(1,N)*2*pi;

Y(1,:) = ro2.*cos(teta2);
Y(2,:) = ro2.*sin(teta2);

figure(3)
hold all
scatter(X(1,:),X(2,:))
scatter(Y(1,:),Y(2,:))
title('Odbirci')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2')

%% Pretpostavka da postoje 3 klase

X = [X Y];

pom = rand(1,2*N);
X1 = [];
X2 = [];
X3 = [];

K = zeros(1,2*N);

for i = 1:2*N
    if pom(i) < 0.33
        K(i)=1;
    elseif pom(i)<0.66
        K(i)=2;
    else 
        K(i)=3;
    end
end

K(1:100) = 1;
K(401:500) = 2;
K(801:900) = 3;

for i=1:2*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
        case 3
            X3 = [X3, X(:,i)];
    end
end

figure()
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
scatter(X3(1,:),X3(2,:));
hold off
title('Slučajna klasterizacija')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2','K3')
%% Klasterizacija sa 3 klase

q = ones(2*N,3)/3;

q(1:100,1) = 1; 
q(401:500,2) = 1;
q(801:900,3) = 1;
q(1:100,2) = 0; 
q(401:500,1) = 0;
q(801:900,1) = 0;
q(801:900,2) = 0;
q(1:100,3) = 0;
q(401:500,3) = 0;

lmax = 1000;
l = 1;

while(l<lmax)
   Pj = zeros(1,3);
   Mj = zeros(2,1,3);
   Sj = zeros(2,2,3);
   fj = zeros(4*N,3);
   
   for j=1:3
       Pj(j) = sum(q(:,j))/(2*N);
 
       Mj(:,:,j)=zeros(2,1);
       for i=1:2*N
           Mj(:,:,j) = Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j) = Mj(:,:,j)/(2*N*Pj(j));
       
 
       Sj(:,:,j)=zeros(2,2);
       for i=1:2*N
           Sj(:,:,j) = Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j) = Sj(:,:,j)/(2*N*Pj(j));
       
   end
   q_pre = q;
   for i=1:2*N
       for j=1:3
           fj(i,j) = 1/(2*pi*det(Sj(:,:,j))^(0.5))*exp(-0.5*(X(:,i)-Mj(:,:,j))'*inv(Sj(:,:,j))*(X(:,i)-Mj(:,:,j)));
       end
       
       for j=1:3
           q(i,j) = Pj(j)*fj(i,j)/(Pj(1)*fj(i,1)+Pj(2)*fj(i,2)+Pj(3)*fj(i,3));    
       end
       
   end
   if max(max(q-q_pre))>10^(-3)
       l = l+1;
   else
       break;
   end
  
end

for i=1:2*N
   [qmax,idx] = max(q_pre(i,:));
   K(i) = idx;
end

figure()
K1 = [];
K2 = [];
K3 = [];

for i=1:2*N
    if(K(i) == 1)
        K1 = [K1, X(:,i)];
    end
    if(K(i) == 2)
        K2 = [K2, X(:,i)];
    end
    if(K(i) == 3)
        K3 = [K3, X(:,i)];
    end
    hold all
end

scatter(K1(1,:),K1(2,:))
hold all;
scatter(K2(1,:),K2(2,:))
scatter(K3(1,:),K3(2,:))
xlabel('x1')
ylabel('x2')
title(['Iteracija broj ', num2str(l)])
legend('K1', 'K2','K3')
%% Generisanje nelinearno separabilnih klasa

close all
clc

N = 500;
ro1 = rand(1,N);
teta1 = rand(1,N)*2*pi;

X = zeros(2,N);
X(1,:) = ro1.*cos(teta1);
X(2,:) = ro1.*sin(teta1);

figure()
hold all
scatter(X(1,:),X(2,:))

ro2 = rand(1,N)+2;
teta2 = rand(1,N)*2*pi;

Y(1,:) = ro2.*cos(teta2);
Y(2,:) = ro2.*sin(teta2);

figure()
scatter(Y(1,:),Y(2,:))
title('Odbirci')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2')

%% Pretpostavka da imamo 4 klase

X = [X Y];

pom = rand(1,2*N);
X1 = [];
X2 =[];
X3 = [];
X4 = [];
K = zeros(1,2*N);

for i = 1:2*N
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
K(251:350) = 2;
K(501:600) = 3;
K(751:850) = 4;

for i=1:2*N
    switch(K(i))
        case 1 
            X1 = [X1, X(:,i)];
        case 2 
            X2 = [X2, X(:,i)];
        case 3 
            X3 = [X3, X(:,i)];
        case 4 
            X4 = [X4, X(:,i)];
    end
end
            
figure()
scatter(X1(1,:), X1(2,:))
hold on
scatter(X2(1,:), X2(2,:))
scatter(X3(1,:), X3(2,:))
scatter(X4(1,:), X4(2,:))
hold off
title('Slučajna klasterizacija')
xlabel('x1')
ylabel('x2')

%% Klasterizacija sa 4 klase

q = ones(2*N,4)/4;

q(1:100,1) = 1; 
q(251:350,2) = 1;
q(501:600,3) = 1;
q(751:850,4) = 1;
q(1:100,2) = 0; 
q(251:350,1) = 0;
q(501:600,1) = 0;
q(751:850,1) = 0;
q(1:100,3) = 0; 
q(251:350,3) = 0;
q(501:600,2) = 0;
q(751:850,2) = 0;
q(1:100,4) = 0; 
q(251:350,4) = 0;
q(501:600,4) = 0;
q(751:850,3) = 0;

lmax = 1000;
l = 1;

while(l<lmax)
   Pj = zeros(1,4);
   Mj = zeros(2,1,4);
   Sj = zeros(2,2,4);
   fj = zeros(2*N,4);
   
   for j=1:4
       Pj(j) = sum(q(:,j))/(2*N);
 
       Mj(:,:,j)=zeros(2,1);
       for i=1:2*N
           Mj(:,:,j) = Mj(:,:,j)+X(:,i)*q(i,j);
       end
       Mj(:,:,j) = Mj(:,:,j)/(2*N*Pj(j));
       
 
       Sj(:,:,j)=zeros(2,2);
       for i=1:2*N
           Sj(:,:,j) = Sj(:,:,j)+(X(:,i)-Mj(:,:,j))*(X(:,i)-Mj(:,:,j))'*q(i,j);     
       end
       Sj(:,:,j) = Sj(:,:,j)/(2*N*Pj(j));
       
   end

   q_pre = q;
   
   for i=1:2*N
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


for i=1:2*N
   [qmax,idx] = max(q_pre(i,:));
   K(i) = idx;
end

figure(5)
K1 = [];
K2 = [];
K3 = [];
K4 = [];

for i=1:2*N
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

