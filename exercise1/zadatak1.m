clear all;
close all;
clc;

kamen_putanja = 'kamen\';
papir_putanja = 'papir\';
makaze_putanja = 'makaze\';

kamen_im_names = dir(fullfile(kamen_putanja, '*.png'));
papir_im_names = dir(fullfile(papir_putanja, '*.png'));
makaze_im_names = dir(fullfile(makaze_putanja, '*.png'));

N1 = length(kamen_im_names);
X1 = zeros(N1, 2);
N2 = length(papir_im_names);
X2 = zeros(N2, 2);
N3 = length(makaze_im_names);
X3 = zeros(N3, 2);

for i=1:N1
    im = imread(fullfile(kamen_putanja, kamen_im_names(i).name));

%     figure();
%     imshow(im);
%     title('Kamen');

    bin = obrada_slike(im);
    X1(i,:) = izdvoji_obelezja(bin);
end
for i=1:N2
    im = imread(fullfile(papir_putanja, papir_im_names(i).name));


%     figure();
%     imshow(im);
%     title('Papir');

    bin = obrada_slike(im);
    X2(i,:) = izdvoji_obelezja(bin);
end
for i=1:N3
    im = imread(fullfile(makaze_putanja, makaze_im_names(i).name));


%     figure();
%     imshow(im);
%     title('Makaze');
    
 
    bin = obrada_slike(im);
    X3(i,:) = izdvoji_obelezja(bin);
end

%%
Ntr = 600;
X1tr = X1(1:Ntr,:);
X1test = X1(Ntr+1:end,:);
X2tr = X2(1:Ntr,:);
X2test = X2(Ntr+1:end,:);
X3tr = X3(1:Ntr,:);
X3test = X3(Ntr+1:end,:);

%%
figure();
plot(X1tr(:,1), X1tr(:,2),'sb', X2tr(:,1), X2tr(:,2),'hg', X3tr(:,1), X3tr(:,2),'ro');
legend('kamen', 'papir', 'makaze');
xlabel('x1');
ylabel('x2');
title('Klase');
grid on;

%% Projektivanje klasifikatora
M1 = mean(X1tr);
M1 = M1';
M2 = mean(X2tr);
M2 = M2';
M3 = mean(X3tr);
M3 = M3';
S1 = cov(X1tr);
S1 = S1';
S2 = cov(X2tr);
S2 = S2';
S3 = cov(X3tr);
S3 = S3';

%% Klasifikacija sve tri klase
Xtest = [X1test; X2test; X3test];
Ylabel = [ones(N1-Ntr, 1); 2*ones(N2-Ntr, 1); 3*ones(N3-Ntr, 1)];
Ypred = ones(length(Ylabel), 1);

Xtest = Xtest';
Ylabel = Ylabel';
Ypred = Ypred';

for i=1:length(Xtest)
    x = Xtest(:,i);
    f1=1/(2*pi*det(S1)^0.5)*exp(-0.5*(x-M1)'*inv(S1)*(x-M1));
    f2=1/(2*pi*det(S2)^0.5)*exp(-0.5*(x-M2)'*inv(S2)*(x-M2));
    f3=1/(2*pi*det(S3)^0.5)*exp(-0.5*(x-M3)'*inv(S3)*(x-M3));
    f = [f1 f2 f3];
    if max(f)==f1
        Ypred(i) = 1;
    else 
        if max(f)==f2
            Ypred(i) = 2;
        else
            Ypred(i) = 3;
        end
    end
end

figure()
cm = confusionchart(Ylabel,Ypred);
val = cm.NormalizedValues;
disp('Tačnost');
disp((val(1,1)+val(2,2)+val(3,3))/sum(sum(val)));

figure();
hold on;
for i=1:length(Xtest)
    if Ypred(i)==1
        if Ylabel(i)==1
            plot(Xtest(1,i), Xtest(2,i), 'sb', 'DisplayName', 'kamen klasifikovan kao kamen');
        end
        if Ylabel(i)==2
            plot(Xtest(1,i), Xtest(2,i), 'hb', 'DisplayName', 'papir klasifikovan kao kamen');
        end
        if Ylabel(i)==3
            plot(Xtest(1,i), Xtest(2,i), 'ob', 'DisplayName', 'makaze klasifikovane kao kamen');
        end
    end
    if Ypred(i)==2
        if Ylabel(i)==1
            plot(Xtest(1,i), Xtest(2,i), 'sg', 'DisplayName', 'kamen klasifikovan kao papir');
        end
        if Ylabel(i)==2
            plot(Xtest(1,i), Xtest(2,i), 'hg', 'DisplayName', 'papir klasifikovan kao papir');
        end
        if Ylabel(i)==3
            plot(Xtest(1,i), Xtest(2,i), 'og', 'DisplayName', 'makaze klasifikovane kao papir');
        end
    end
    if Ypred(i)==3
        if Ylabel(i)==1
            plot(Xtest(1,i), Xtest(2,i), 'sr', 'DisplayName', 'kamen klasifikovan kao makaze');
        end
        if Ylabel(i)==2
            plot(Xtest(1,i), Xtest(2,i), 'hr', 'DisplayName', 'papir klasifikovan kao makaze');
        end
        if Ylabel(i)==3
            plot(Xtest(1,i), Xtest(2,i), 'or', 'DisplayName', 'makaze klasifikovane kao makaze');
        end
    end
end
title('Klasifikacije');
grid on;

%% Histogrami
Xtest = [X1test; X3test];
Ylabel = [ones(N1-Ntr, 1); 2*ones(N3-Ntr, 1)];
Ypred = ones(length(Ylabel), 1);
Xtest = Xtest';
Ylabel = Ylabel';
Ypred = Ypred';

figure();
subplot(2,1,1);
histogram(Xtest(1,:));
title('Histogram prvog obeležja');
subplot(2,1,2);
histogram(Xtest(2,:));
title('Histogram drugog obeležja');

%% Linearni klasifikator na bazi zeljenog izlaza za kamen i makaze
Xtest = [X1test; X3test];
Ylabel = [ones(N1-Ntr, 1); 2*ones(N3-Ntr, 1)];
Ypred = ones(length(Ylabel), 1);
Xtest = Xtest';
Ylabel = Ylabel';
Ypred = Ypred';

U = [-1*ones(1,Ntr) ones(1,Ntr); -X1tr', X3tr'];
G = [ones(Ntr,1);  1.7*ones(Ntr,1)]; 
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3);

v_opt = [v1; v2];
for i =1 : length(Xtest)
   x = Xtest(:,i);
   if (v_opt' * x + v0) > 0
       Ypred(i) = 2;  
   end  
end

x1 = 0:1:200;
x2 = -(v0+v1*x1)/v2;

figure();
hold on;
for i=1:length(Xtest)
    if Ypred(i)==1
        if Ylabel(i)==1
            plot(Xtest(1,i), Xtest(2,i), 'sb', 'DisplayName', 'kamen klasifikovan kao kamen');
        end
        if Ylabel(i)==2
            plot(Xtest(1,i), Xtest(2,i), 'hb', 'DisplayName', 'makaze klasifikovane kao kamen');
        end
    end
    if Ypred(i)==2
        if Ylabel(i)==1
            plot(Xtest(1,i), Xtest(2,i), 'sg', 'DisplayName', 'kamen klasifikovan kao makaze');
        end
        if Ylabel(i)==2
            plot(Xtest(1,i), Xtest(2,i), 'hg', 'DisplayName', 'makaze klasifikovane kao makaze');
        end
    end
end
plot(x1,x2,'k--');
title('Klasifikacije kamena i makaza');
grid on;


C = confusionmat(Ylabel,Ypred);
figure()
cm = confusionchart(Ylabel,Ypred);
val = cm.NormalizedValues;
disp('Tačnost');
disp((val(1,1)+val(2,2))/sum(sum(val)));