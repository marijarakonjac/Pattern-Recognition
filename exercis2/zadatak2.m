clear
close all
clc

%% Kreiranje i prikaz odbiraka klasa

N = 500;

M11 = [2, 1]';
S11 = [0.7 0.3; 0.3 0.7];
M12 = [1, 3]';
S12 = [0.9 0.5; 0.5 0.9];

M21 = [-2, 0]';
S21 = [0.8 0.3; 0.3 0.8];
M22 = [-2.5, 2.5]';
S22 = [1 0.5; 0.5 1];

P11 = 0.6;
P12 = 0.4;
P21 = 0.7;
P22 = 0.3;

K11 = mvnrnd(M11, S11, N)';
K12 = mvnrnd(M12, S12, N)';
K21 = mvnrnd(M21, S21, N)';
K22 = mvnrnd(M22, S22, N)';

pom1 = rand(1, N);
pom2 = rand(1, N);

K1 = (pom1 < P11).*K11 + (pom1 >= P11).*K12;
K2 = (pom2 < P21).*K21 + (pom2 >= P21).*K22;

figure(1)
scatter(K1(1, :), K1(2, :), 'ro')
hold all
scatter(K2(1, :), K2(2, :), 'bo')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2')
title('Odbirci klasa')
hold off

%% Teorijska procena fgv

x = -8:0.1:6;
y = -5:0.1:9;

f1 = zeros(length(x),length(y));
f2 = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i); y(j)];
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
        f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
        f1(i,j) = P11*f11 + P12*f12;
        f2(i,j) = P21*f21 + P22*f22; 
    end
end

figure(2)
mesh(x, y, f1')
hold all
mesh(x, y, f2')
hold off
xlabel('x1');
ylabel('x2');
title('Teorijska funkcija gustine verovatnoće');

%% Histogrami

figure(3)
bins = 20; 
hist3([K1'; K2'],'Nbins', [bins, bins],'CDataMode', 'auto', 'FaceColor', 'interp')
xlim([-8 6])
ylim([-5 9])
xlabel('x1')
ylabel('x2')
title('Histogram generisanih odbiraka')


%% Bajesov klasifikator minimalne greške

h = -log(f1./f2);

figure(4)
scatter(K1(1, :), K1(2, :), 'ro')
hold all
scatter(K2(1, :), K2(2, :), 'bo')		
contour(x, y, h',[0 0], 'k', 'LineWidth', 1.2)
title('Klasifikacija Bajesovim klasifikatorom')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2', 'Diskriminaciona kriva')

% Poređenje grešaka 

Xu = [K1, K2];
Y_true = [ones(1, N), 2*ones(1, N)];
Y_pred = zeros(1, 2*N);

for i = 1:length(Xu)
    Xn = Xu(:, i);
    f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(Xn-M11)'*inv(S11)*(Xn-M11));
    f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(Xn-M12)'*inv(S12)*(Xn-M12));
    f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(Xn-M21)'*inv(S21)*(Xn-M21));    
    f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(Xn-M22)'*inv(S22)*(Xn-M22));
    f1n = P11*f11 + P12*f12;
    f2n = P21*f21 + P22*f22; 
    
    if(f1n > f2n)
        Y_pred(i) = 1;
    else
        Y_pred(i) = 2;
    end
end

C = confusionmat(Y_true, Y_pred);
disp(C)
e1n = C(1,2)/N;
e2n = C(2,1)/N;

%Teorijska procena greške

e1t = 0;
e2t = 0;

for i = 1:length(x)
    for j = 1:length(y)
        if(h(i, j) > 0) 
            e1t = e1t + f1(i,j)*0.1*0.1;   
        else
            e2t = e2t + f2(i,j)*0.1*0.1;  
        end
    end
end

disp('Greška prvog tipa (eksperimentalno i teorijski): ')
s1 = sprintf('\t \t %.3f \t %.3f \n', e1n, e1t);
disp(s1)
disp('Greška drugog tipa (eksperimentalno i teorijski): ')
s2 = sprintf('\t \t %.3f \t %.3f \n', e2n, e2t);
disp(s2)

%% Klasifikator minimalne cene

c11 = 0;
c22 = 0;
c12 = 1; 
c21 = 5;

k = -log(f1./f2);  
t = -log((c12 - c22) / (c21 - c11));

contour(x, y, k',[t t], 'g--', 'LineWidth', 1.2)
legend('K1', 'K2', 'Bajesov klasifikator', 'Klasifikator minimalne cene')
hold off

%% Neuman-Pearsonov klasifikator

e0 = e1n + e2n;

%različite vrednosti parametra mi

br = 0;
for mi = 0.01:0.01:10
    br = br + 1;
    e2np(br) = 0;
    for i = 1:length(x)
        for j = 1:length(y)
            if (h(i,j) <= -log(mi))
                e2np(br) = e2np(br) + 0.1*0.1*f2(i,j); 
            end
        end
    end
end

figure(5)
plot(0.01:0.01:10, e2np, 'k')
hold all
yline(e0, '-.r')
xlabel('\mu')
ylabel('\epsilon_2')
title('Zavisnost greške drugog tipa od \mu')
legend('Zavisnost greške drugog tipa od \mu', '\epsilon_0')

%Klasifikacija Neuman-Pearsonovim klasifikatorom i trazenje mi za zadato e0

[val, idx] = min(abs(e2np - e0));
miv = e2np(idx);
t = -log(miv);

figure(6)
scatter(K1(1, :), K1(2, :), 'ro')
hold all
scatter(K2(1, :), K2(2, :), 'bo')		
contour(x, y, h',[t t], 'k', 'LineWidth', 1.2)
title('Klasifikacija NP klasifikatorom')
xlabel('x1')
ylabel('x2')
legend('K1', 'K2', 'Klasifikaciona linija');

%% Wald-ov sekvencijalni test

e1 = 10^-10;
e2 = 10^-10;
a = -log((1-e1)/e2);
b = -log(e1/(1-e2));

Xp = [K1, K2];

for i = 1:100
    Sm = 0;
    Sm_arr = [];
    br = 1;
    while (Sm > a && Sm < b)
        idx = randi(size(Xp, 2));
        X = Xp(:, idx);
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
        f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
        f1 = P11*f11 + P12*f12;
        f2 = P21*f21 + P22*f22;
        h = -log(f1/f2);
        Sm = Sm + h;
        Sm_arr(br) = Sm;
        br = br + 1;
        
        if Sm <= a
            figure(7)
            plot(Sm_arr(1:br-1), 'r')
            hold all
        elseif Sm >= b
            figure(7)
            plot(Sm_arr(1:br-1), 'g')
            hold all
        end
    end
end

figure(7)
yline(a, 'k-')
yline(b, 'k-')
xlabel('m')
ylabel('s(m)')
title('Wald-ov sekvencijalni test')

% Zavisnost broja potrebnih odbiraka prve klase od epsilon 1 i 2

e1_arr = logspace(-10, 0, 100);
e2_arr = logspace(-10, 0, 100);
m1_arr = zeros(1,length(e1_arr));

for i = 1:length(e1_arr)
    m1_e = [];
    a1 = -log((1-e1_arr(i))/e2);
    b1 = -log(e1_arr(i)/(1-e2));

    for j = 1:100
        Sm = 0;
        Sm_arr = [];
        br = 1;

        while (Sm > a1 && Sm < b1)
           idx = randi(size(K1,2));
           X = K1(:,idx);
           f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
           f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
           f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
           f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
           f1 = P11*f11 + P12*f12;
           f2 = P21*f21 + P22*f22;
           h = -log(f1/f2);
           Sm = Sm + h;
           Sm_arr(br) = Sm;
           br = br + 1;

           if Sm > a1
              m1_e(end+1) = br-1;
           end

        end
    
    end
    m1_arr(i) = mean(m1_e);
end

m2_arr = zeros(1,length(e2_arr));

for i = 1:length(e2_arr)
    m2_e = [];
    a2 = log(e2_arr(i)/(1-e1));
    b2 = log((1-e2_arr(i))/e1);
    for j = 1:100
        Sm = 0;
        Sm_arr = [];
        br = 1;

        while (Sm > a2 && Sm < b2)
           idx = randi(size(K1,2));
           X = K1(:,idx);
           f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
           f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
           f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
           f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
           f1 = P11*f11 + P12*f12;
           f2 = P21*f21 + P22*f22;
           h = -log(f1/f2);
           Sm = Sm + h;
           Sm_arr(br) = Sm;
           br = br + 1;

           if Sm <= a2
              m2_e(end+1) = br-1;
           end

        end
    
    end
    m2_arr(i) = mean(m2_e);
end

figure(8)
plot(e1_arr, m1_arr)
set(gca, 'XScale', 'log')
hold all
plot(e2_arr, m2_arr)
set(gca, 'XScale', 'log')
xlim([10^-10, 1]) 
xlabel('\epsilon_1 ili \epsilon_2')
ylabel('m')
title('Zavisnost broja potrebnih odbiraka iz prve klase od \epsilon_1 i \epsilon_2')
legend('Zavisnost od \epsilon_1', 'Zavisnost od \epsilon_2')

% Zavisnost broja potrebnih odbiraka druge klase od epsilon 1 i 2
m1_arr = zeros(1,length(e1_arr));

for i = 1:length(e1_arr)
    m1_e = [];
    a1 = -log((1-e1_arr(i))/e2);
    b1 = -log(e1_arr(i)/(1-e2));
    for j = 1:100
        Sm = 0;
        Sm_arr = [];
        br = 1;

        while (Sm > a1 && Sm < b1)
           idx = randi(size(K2,2));
           X = K2(:,idx);
           f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
           f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
           f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
           f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
           f1 = P11*f11 + P12*f12;
           f2 = P21*f21 + P22*f22;
           h = -log(f1/f2);
           Sm = Sm + h;
           Sm_arr(br) = Sm;
           br = br + 1;

           if Sm >= b1
              m1_e(end+1) = br-1;
           end

        end
    
    end
    m1_arr(i) = mean(m1_e);
end

m2_arr = zeros(1,length(e2_arr));

for i = 1:length(e2_arr)
    m2_e = [];
    a2 = log(e2_arr(i)/(1-e1));
    b2 = log((1-e2_arr(i))/e1);
    for j = 1:100
        Sm = 0;
        Sm_arr = [];
        br = 1;

        while (Sm > a2 && Sm < b2)
           idx = randi(size(K2,2));
           X = K2(:,idx);
           f11 = 1/(2*pi*det(S11)^0.5)*exp(-1/2*(X-M11)'*inv(S11)*(X-M11));
           f12 = 1/(2*pi*det(S12)^0.5)*exp(-1/2*(X-M12)'*inv(S12)*(X-M12));
           f21 = 1/(2*pi*det(S21)^0.5)*exp(-1/2*(X-M21)'*inv(S21)*(X-M21));    
           f22 = 1/(2*pi*det(S22)^0.5)*exp(-1/2*(X-M22)'*inv(S22)*(X-M22));
           f1 = P11*f11 + P12*f12;
           f2 = P21*f21 + P22*f22;
           h = -log(f1/f2);
           Sm = Sm + h;
           Sm_arr(br) = Sm;
           br = br + 1;

           if Sm >= b2
              m2_e(end+1) = br-1;
           end

        end
    
    end
    m2_arr(i) = mean(m2_e);
end

figure(9)
plot(e1_arr, m1_arr)
set(gca, 'XScale', 'log')
hold all
plot(e2_arr, m2_arr)
set(gca, 'XScale', 'log')
xlim([10^-10, 1]) 
xlabel('\epsilon_1 ili \epsilon_2')
ylabel('m')
title('Zavisnost broja potrebnih odbiraka iz druge klase od \epsilon_1 i \epsilon_2')
legend('Zavisnost od \epsilon_1', 'Zavisnost od \epsilon_2')

