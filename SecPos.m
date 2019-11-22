Ke = 0:100;
Km = 32;
m = (Km+1).^2 ./ (2.*Km+1);
n = (Ke+1).^2 ./ (2.*Ke+1);
SIGMA = 2.*Km+1;
PHI = 2.*Ke+1;
z =(n.*SIGMA) ./ (m.*PHI) ;
y = (n.*beta(m,n)).^-1;
x =  hyp2f1(n,m+n,n+1, -z);


w =  hyp2f1(n,m+n,n+1, -ones(1,size(z,2)));

v =  n.*beta(m,n);

Ph = z.^n .* y .* x

plot (Ke/Km,w)
figure;plot (Ke/Km,v)

%plot(Ke/Km,Ph)
%figure;plot (Ke/Km,x)
%figure;plot (Ke/Km,z)
La_s = 1e-3:1e-3:1e-1;
La_e = 60e-6;
Pr = La_s ./ (La_s + La_e);


%figure;plot(La_s,Pr)
hold on
La_e = 60e-6;

Pr = La_s ./ (La_s + La_e);


%plot(La_s,Pr)