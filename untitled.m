    dd =1e-4;
    do =1e-4;
    dw=10;
    b0=1./(1+exp(do/dw));
    b1=1./(1+exp((do-dd)/dw));
    dd=(b1-b0)*dw*4;
    [b0 b1],dd
