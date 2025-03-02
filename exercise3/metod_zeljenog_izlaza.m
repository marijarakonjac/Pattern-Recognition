function [v0,V] = metod_zeljenog_izlaza(X1,X2,gama)
N1 = length(X1(1,:));
N2 = length(X2(1,:));

Z = [-1*ones(1,N1) 1*ones(1,N2);-X1 X2];
W =(Z*Z')^(-1)*Z * gama;
v0 = W(1); 
V = W(2:3);

end