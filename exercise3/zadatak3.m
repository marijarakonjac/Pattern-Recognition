clear
close all
clc

%% Generisanje klasa i prikaz klasa

M1 = [0 4]'; S1 = [0.9 0.8;0.8 0.9];
M2 = [8 4]'; S2 = [0.8 -0.5;-0.5 0.8];
M3 = [4 -1]'; S3 = [0.5 0.3;0.3 0.5];

N = 1000;

X1 =  mvnrnd(M1,S1,N)';
X2 =  mvnrnd(M2,S2,N)';
X3 =  mvnrnd(M3,S3,N)';

figure()
plot(X1(1,:),X1(2,:),'ro',X2(1,:),X2(2,:),'bx',X3(1,:),X3(2,:),'g*')
title('Generisane klase');
grid on;
xlabel('x_1'); 
ylabel('x_2');
legend('K1','K2','K3'); 

%% Projektovanje linearnog klasifikatora metodom resupstitucije

X1_novo = [X1,X2];
X2_novo = X3;
X_test = [X1 X2 X3];

[s_opt, v0_opt, Neps_opt, M1_est, M2_est, S1_est, S2_est] = resupstitucija(X1_novo, X2_novo);

V = (s_opt * S1_est + (1-s_opt) * S2_est)^(-1) * (M2_est - M1_est);
v0 = v0_opt;

x1 = -5:0.1:10;
x2 = -(v0 + V(1) * x1) / V(2);
klasif1 = [x1; x2];

Y = V' * X_test + v0;
Y(Y > 0) = 3;
indeksi = find(Y < 0);

X1_novo = X1;
X2_novo = X2;

[s_opt, v0_opt, Neps_opt, M1_est, M2_est, S1_est, S2_est] = resupstitucija(X1_novo, X2_novo);

V = (s_opt * S1_est + (1-s_opt) * S2_est)^(-1) * (M2_est - M1_est);
v0 = v0_opt;

x1 = -5:0.1:10;
x2 = -(v0 + V(1) * x1) / V(2);
klasif2 = [x1; x2];

if ~isempty(indeksi)
    Y(indeksi) = V' * X_test(:,indeksi) + v0;
    Y(Y < 0) = 1;
    Y((Y ~= 1) & (Y ~= 3)) = 2;
end

figure()
plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bx', X3(1,:), X3(2,:), 'g*', klasif1(1,:), klasif1(2,:), klasif2(1,:), klasif2(2,:))
legend('K1', 'K2', 'K3', 'Klasifikator 12/3', 'Klasifikator 1/2');
title('Prikaz odbiraka i klasifikatora');
xlabel('x_1'); 
ylabel('x_2');
ylim([-4 12]);
grid on;

Y_prava = [ones(1, N) 2 * ones(1, N) 3 * ones(1, N)];
c = confusionmat(Y_prava, Y);
disp(c)

figure()
confusionchart(c)

%% Projektovanje linearnog klasifikatora metodom zeljenog izlaza

%gama1 = [2 * ones(2 * N, 1); -ones(N, 1)];
%gama = [ones(3 * N, 1) 2 * ones(3 * N, 1) gama1]';

tezine = [1 2 3];

for i = 1:3 
    gama1 = [ones(2 * N, 1); tezine(i)*ones(N, 1)];
    gama2 = [ones(N, 1); tezine(i)*ones(N, 1)];

    X1_novo = [X1, X2];
    X2_novo = X3;

    [v0, V] = metod_zeljenog_izlaza(X1_novo, X2_novo, gama1);

    x1 = -5:0.1:10;
    x2 = (-v0 - V(1) * x1) / V(2);
    klasif1 = [x1; x2];
    
    Y = V' * X_test + v0;
    Y(Y > 0) = 3;
    indeksi = find(Y < 0);

    X1_novo = X1;
    X2_novo = X2;

    [v0, V] = metod_zeljenog_izlaza(X1_novo, X2_novo, gama2);

    x1 = -5:0.1:10;
    x2 = (-v0 - V(1) * x1) / V(2);
    klasif2 = [x1; x2];
    
    if ~isempty(indeksi)
        Y(indeksi) = V' * X_test(:,indeksi) + v0;
        Y(Y < 0) = 1;
        Y((Y ~= 1) & (Y ~= 3)) = 2;
    end

    figure()
    plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bx', X3(1,:), X3(2,:), 'g*', klasif1(1,:), klasif1(2,:), klasif2(1,:), klasif2(2,:))
    legend('K1', 'K2', 'K3', 'Klasifikator 12/3', 'Klasifikator 1/2');
    xlabel('x_1'); 
    ylabel('x_2');
    ylim([-4 12]);
    title(['Prikaz odbiraka i klasifikatora za gama ', num2str(i)]);
    grid on;
            
    Y_prava = [ones(1, N) 2 * ones(1, N) 3 * ones(1, N)];
    c = confusionmat(Y_prava, Y);
    disp(c)

    figure()
    confusionchart(c)
    
end

%% Projektovanje kvadratnog klasifikatora

N = 1000;

M1 = [0;0];
R1 = 2 * rand(1,N);
u1 = 2 * pi * rand(1,N);
X1 = [R1.*cos(u1) ; R1.*sin(u1)] + M1 * ones(1, N);

M2 = [0;0];
Ru = 2.5;
Rd = 1.5;
R2 = Rd * rand(1, N) + Ru;
u2 = 2 * pi * rand(1, N);
X2 = [R2.*cos(u2) ; R2.*sin(u2)] + M2 * ones(1, N);

figure
plot(X1(1,:), X1(2,:), 'r*', X2(1,:), X2(2,:), 'b*')
legend('Klasa 1', 'Klasa 2'); 
xlabel('x_1'); 
ylabel('x_2');
title('Nelinearno separabilne klase separabilnih klasa');

gama = ones(2 * N, 1);

Z = [-1 * ones(1, N) 1 * ones(1, N); -X1 X2; -X1(1,:).^2 X2(1,:).^2; -2 * X1(1,:).*X1(2,:) 2 * X2(1,:).*X2(2,:); -X1(2,:).^2 X2(2,:).^2];
W = (Z * Z')^(-1) * Z * gama;
v0 = W(1); V = W(2:3);
Q = [W(4) W(5); W(5) W(6)];

syms xp yp 

[xp, yp, ~, ~] = solve(v0 + xp * V(1) + yp * V(2) + xp^2 * Q(1) + 2 * xp * yp * Q(2) + yp^2 * Q(4), xp, yp, 'returnCondition', true);
z = -3:0.001:5;
xp = eval(xp);
xp = [xp(1,:) fliplr(xp(2,:))];
xp1 = xp(imag(xp) == 0);
yp = eval(yp);
yp = [yp(1,:) fliplr(yp(2,:))];
yp1 = yp(imag(xp) == 0);

figure()
hold on;
plot(X1(1,:), X1(2,:), 'r*', X2(1,:), X2(2,:), 'b*')
plot(xp1, yp1, 'k')
legend('Klasa 1', 'Klasa 2', 'Klasifikator');
xlabel('x_1'); 
ylabel('x_2');
title('Prikaz odbiraka i diskriminacione krive');
