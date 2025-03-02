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

figure()
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

%% C mean klasterizacija

X = [K1' K2' K3' K4'];

%inicijalna klasterizacija
pom = rand(1,4*N);

X1 = [];
X2 = [];
X3 = [];
X4 = [];

for i = 1:4*N
    if pom(i) < 0.25
        X1 = [X1, X(:, i)];
    elseif pom(i) < 0.5
        X2 = [X2, X(:, i)];
    elseif pom(i) < 0.75
        X3 = [X3, X(:,i)];
    else 
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

%% Primena algoritma klasterizacije

lmax = 100;  
reklas = 1;  

M1 = mean(X1, 2);
M2 = mean(X2, 2);
M3 = mean(X3, 2);
M4 = mean(X4, 2);

N1 = length(X1(1,:));
N2 = length(X2(1,:));
N3 = length(X3(1,:));
N4 = length(X4(1,:));

l = 1;

while(l < lmax) && reklas

    X1pom = [];
    X2pom = [];
    X3pom = [];
    X4pom = [];

    reklas = 0;

    for i = 1:N1
        d1 = (X1(1,i)-M1(1,1))^2+(X1(2,i)-M1(2,1))^2;
        d2 = (X1(1,i)-M2(1,1))^2+(X1(2,i)-M2(2,1))^2;
        d3 = (X1(1,i)-M3(1,1))^2+(X1(2,i)-M3(2,1))^2;
        d4 = (X1(1,i)-M4(1,1))^2+(X1(2,i)-M4(2,1))^2;
        
        dmin = min([d1,d2,d3,d4]);
        if dmin == d1
            X1pom = [X1pom, X1(:, i)];
        elseif dmin == d2
            X2pom = [X2pom, X1(:, i)];  
            reklas = 1;
        elseif dmin == d3
            reklas = 1;
            X3pom = [X3pom, X1(:,i)];
        else 
            reklas = 1;
            X4pom = [X4pom, X1(:,i)];
        end
    end
    for i = 1:N2
        d1 = (X2(1,i)-M1(1,1))^2+(X2(2,i)-M1(2,1))^2;
        d2 = (X2(1,i)-M2(1,1))^2+(X2(2,i)-M2(2,1))^2;
        d3 = (X2(1,i)-M3(1,1))^2+(X2(2,i)-M3(2,1))^2;
        d4 = (X2(1,i)-M4(1,1))^2+(X2(2,i)-M4(2,1))^2;
        
        dmin = min([d1,d2,d3,d4]);
        if dmin == d1
            X1pom = [X1pom, X2(:, i)];
            reklas = 1;
        elseif dmin == d2
            X2pom = [X2pom, X2(:, i)];  
        elseif dmin == d3
            reklas = 1;
            X3pom = [X3pom, X2(:,i)];
        else 
            reklas = 1;
            X4pom = [X4pom, X2(:,i)];
        end
    end
    for i=1:N3
        d1 = (X3(1,i)-M1(1,1))^2+(X3(2,i)-M1(2,1))^2;
        d2 = (X3(1,i)-M2(1,1))^2+(X3(2,i)-M2(2,1))^2;
        d3 = (X3(1,i)-M3(1,1))^2+(X3(2,i)-M3(2,1))^2;
        d4 = (X3(1,i)-M4(1,1))^2+(X3(2,i)-M4(2,1))^2;
        
        dmin = min([d1,d2,d3,d4]);
        if dmin == d1
            X1pom = [X1pom, X3(:, i)];
            reklas = 1;
        elseif dmin == d2
            X2pom = [X2pom, X3(:, i)]; 
            reklas = 1;
        elseif dmin == d3
            X3pom = [X3pom, X3(:,i)];
        else 
            reklas = 1;
            X4pom = [X4pom, X3(:,i)];
        end
    end
    for i=1:N4
        d1 = (X4(1,i)-M1(1,1))^2+(X4(2,i)-M1(2,1))^2;
        d2 = (X4(1,i)-M2(1,1))^2+(X4(2,i)-M2(2,1))^2;
        d3 = (X4(1,i)-M3(1,1))^2+(X4(2,i)-M3(2,1))^2;
        d4 = (X4(1,i)-M4(1,1))^2+(X4(2,i)-M4(2,1))^2;
        
        dmin= min([d1,d2,d3,d4]);
        if dmin==d1
            X1pom = [X1pom, X4(:, i)];
            reklas=1;
        elseif dmin==d2
            X2pom = [X2pom, X4(:, i)]; 
            reklas = 1;
        elseif dmin==d3
            reklas=1;
            X3pom=[X3pom, X4(:,i)];
        else 
            X4pom=[X4pom, X4(:,i)];
        end
    end

    clear X1 X2 X3 X4

    X1 = X1pom;
    X2 = X2pom;
    X3 = X3pom;
    X4 = X4pom;
    
    if(~isempty(X1))
    N1 = length(X1(1,:));
    M1 = mean(X1, 2);
    else
        N1 = 0;
        M1 = [0;0];
    end
    if(~isempty(X2))
    N2 = length(X2(1,:));
    M2 = mean(X2, 2);
    else 
        N2=0; M2=[0;0];
    end
    if(~isempty(X3))
    N3 = length(X3(1,:));
    M3 = mean(X3,2);
    else
        N3 = 0;
        M3 = [0;0];
    end
    if(~isempty(X4))
        N4 = length(X4(1,:));
        M4 = mean(X4,2);
    else
        N4 = 0;
        M4 = [0;0];
    end
    l = l+1;
end

figure()
scatter(X1(1,:), X1(2,:)) 
hold on
scatter(X2(1,:), X2(2,:))
scatter(X3(1,:), X3(2,:))
scatter(X4(1,:), X4(2,:))
title(['Iteracija broj ' num2str(l)])
legend('K1','K2','K3','K4');
xlabel('x1')
ylabel('x2')
