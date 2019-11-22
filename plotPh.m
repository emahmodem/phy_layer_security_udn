Ke = 0:100;
Km = 32;
m = (Km+1).^2 ./ (2.*Km+1);
n = (Ke+1).^2 ./ (2.*Ke+1);
SIGMA = 2.*Km+1;
PHI = 2.*Ke+1;

a = [1/2 1 2 5 10];
Ph = zeros(numel(a),numel(Ke));


for i = 1:numel(a)    
    z =(n.*SIGMA) ./ (a(i) * m.*PHI) ;
    y = (n.*beta(m,n)).^-1;
    x =  hyp2f1(n,m+n,n+1, -z);
    Ph(i,:) = z.^n .* y .* x;
end
%w =  hyp2f1(n,m+n,n+1, -ones(1,size(z,2)));

%v =  n.*beta(m,n);



plot (Ke/Km,Ph)

