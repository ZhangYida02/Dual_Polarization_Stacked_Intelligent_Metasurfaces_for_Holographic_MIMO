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

GTilde00 = sqrt((1-epsilon)*rho_sq/2) * (randn(N1*N2,M1*M2) + 1j*randn(N1*N2,M1*M2)); 
GTilde10 = sqrt(epsilon*rho_sq/2) * (randn(N1*N2,M1*M2) + 1j*randn(N1*N2,M1*M2)); 
GTilde11 = sqrt((1-epsilon)*rho_sq/2) * (randn(N1*N2,M1*M2) + 1j*randn(N1*N2,M1*M2)); 
GTilde01 = sqrt(epsilon*rho_sq/2) * (randn(N1*N2,M1*M2) + 1j*randn(N1*N2,M1*M2)); 

G00=sqrtm(RRx)*GTilde00*sqrtm(RTx);
G10=sqrtm(RRx)*GTilde10*sqrtm(RTx);
G11=sqrtm(RRx)*GTilde11*sqrtm(RTx);
G01=sqrtm(RRx)*GTilde01*sqrtm(RTx);

DpG=[G00 G01;G10 G11];
%%
%Random
% PL=20*log10(4*pi/lambda)+10*3.5*log10(zRx)+9;
% rho_sq=1/db2pow(PL);
% DpG = sqrt(rho_sq/2) * (randn(2*N1*N2,2*M1*M2) + 1j*randn(2*N1*N2,2*M1*M2)); 