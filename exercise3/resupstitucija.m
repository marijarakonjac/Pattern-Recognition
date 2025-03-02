function [s_opt,v0_opt,Neps_opt,M1_est,M2_est,S1_est,S2_est] = resupstitucija(X1,X2)

N1 = length(X1(1,:));
N2 = length(X2(1,:));

M1_est = sum(X1')/N1;
M1_est = M1_est' ;
M2_est = sum(X2')/N2;
M2_est = M2_est' ;

S1_est = zeros(2,2);
S2_est = zeros(2,2);

for i=1:N1
    S1_est = S1_est + (X1(:,i)- M1_est)*(X1(:,i)- M1_est)' ;
end

for i=1:N2
    S2_est = S2_est + (X2(:,i)- M2_est)*(X2(:,i)- M2_est)' ;
end

S1_est = S1_est/N1;
S2_est = S2_est/N2;

s = 0:0.001:1;
v0_opt_s = [];
Neps_s = [];

for i = 1 : length(s)
    V = (s(i)*S1_est+ (1-s(i))*S2_est )^(-1) * (M2_est - M1_est);
    Y1 = V' * X1;
    Y2 = V' * X2;
    Y = [Y1,Y2];
    Y = sort(Y);
    v0 = [];
    Neps = [];
    for j = 1:length(Y)-1
        v0(j) = -(Y(j)+Y(j+1))/2;
        Neps(j)=0;
        for k = 1:N1
            if Y1(k) > -v0(j)
                Neps(j)= Neps(j)+1;
            end         
        end
        for k = 1:N2
            if Y2(k) < -v0(j)
                Neps(j)= Neps(j)+1;
            end         
        end  
        
    end 
    Neps_s(i) = min(Neps);
    v0_opt_s(i) = v0(find(Neps == Neps_s(i),1));     
end

Neps_opt = min(Neps_s);
v0_opt = v0_opt_s(find(Neps_s == Neps_opt,1));
s_opt = s(find(Neps_s == Neps_opt,1));

end