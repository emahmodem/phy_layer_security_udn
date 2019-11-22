

La_s = 1e-3:1e-3:1e-1;
La_e = 60e-6;

Pr = zeros(4,numel(La_s));

alpha = 4;
a = 1;
Pr(1,:) = La_s ./ (La_s + a^(2/alpha)*La_e);

alpha = 4;
a = 2;
Pr(2,:) = La_s ./ (La_s + a^(2/alpha)*La_e);


alpha = 4;
a = 5;
Pr(3,:) = La_s ./ (La_s + a^(2/alpha)*La_e);

alpha = 4;
a = 10;
Pr(4,:) = La_s ./ (La_s + a^(2/alpha)*La_e);

plot(La_s,Pr)
