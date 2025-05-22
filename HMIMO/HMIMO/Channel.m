%%channel
RTx=zeros(M1*M2, M1*M2);
RRx=zeros(N1*N2, N1*N2);
GTilde=zeros(N1*N2, M1*M2);

for m1= 1:M1*M2
    for m2 = 1:M1*M2
        [rowl1,coll1]=ind2sub([M1 M2],m1);
        [rowl2,coll2]=ind2sub([M1 M2],m2);
        r(m1,m2)=sqrt((XTxSim(rowl1)-XTxSim(rowl2))^2+(YTxSim(coll1)-YTxSim(coll2))^2);
        RTx(m1,m2) = sinc(2*r(m1,m2)/lambda);
    end
end

for n1 = 1:N1*N2
    for n2 = 1:N1*N2
        [rowl1,coll1]=ind2sub([N1 N2],n1);
        [rowl2,coll2]=ind2sub([N1 N2],n2);
        t(n1,n2)=sqrt((XRxSim(rowl1)-XRxSim(rowl2))^2+(YRxSim(coll1)-YRxSim(coll2))^2);
        RRx(n1,n2) = sinc(2*t(n1,n2)/lambda);
    end
end

PL=20*log10(4*pi/lambda)+10*3.5*log10(zRx)+9;
rho_sq=1/db2pow(PL);

GTilde = sqrt((1-alp)*rho_sq/2) * (randn(N1*N2,M1*M2) + 1j*randn(N1*N2,M1*M2)); 

G=sqrtm(RRx)*GTilde*sqrtm(RTx);

