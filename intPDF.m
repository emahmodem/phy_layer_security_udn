function intPDF
    x = 0.0000005:0.000001:0.001;
    a = 4;
    LA = 0.001;
    fI = zeros(1,numel(x));
    fI_approx = zeros(1,numel(x));
    fI_approx_2 = zeros(1,numel(x));
    
    for i = 1:numel(x)
        fI(i) = getPDF(x(i),a,LA);  
    end
    fI_approx = pi * LA * (2/a) * gamma(1 + (2/a)) * x.^-(1 + (2/a));
    fI_approx_2 = pi * LA * (2/a) * gamma(1 + (2/a)) * x.^-(1 + (4/a));
    %plot(x,fI,'r-',x,fI_approx,'b.',x,fI_approx_2,'k--')
    %plot(x,fI,'r-',x,fI_approx_2,'k--')
    plot(x,fI,'r-',x,fI_approx,'b-.')
    figure; plot(x,fI,'r-')
end






function f_I = getPDF(x,a,LA)
    k = 1:100;
    A = (-1).^(k+1);
    B = gamma(1 + (2/a) * k);
    C = sin((2/a) * pi * k);
    D = factorial(k);

    E = (pi * LA * gamma(1 + (2/a)) * gamma(1 - (2/a)) / x^(2/a)).^k;

    f_I = 1 / (pi * x) * sum(A.*B.*C.*E./D);
end