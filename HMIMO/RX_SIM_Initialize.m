
%%TX_SIM
XRxSim=-lambda/2*(N1-1)/2:lambda/2:lambda/2*(N1-1)/2;
YRxSim=-lambda/2*(N2-1)/2:lambda/2:lambda/2*(N2-1)/2;
ZRxSim=zRx:-dzRx:zRx-dzRx*(K-1);
WRxSim=zeros(N1*N2,N1*N2,K-1);
DRxSim=zeros(N1*N2,N1*N2,K-1);
ThetaRxSim=zeros(N1*N2,N1*N2,K-1);


%%RxStream
DRxStream=zeros(S,N1*N2);
WRxStream=zeros(S,N1*N2);
ThetaRxStream=zeros(S,N1*N2);
XRxStream=-lambda/2*(S-1)/2:lambda/2:lambda/2*(S-1)/2;
YRxStream=zeros(1,S);
ZRxStream=zRx+dzRx*ones(1,S);
[X_gird,Y_gird]=meshgrid(XRxStream,YRxStream);


%%TX_SIM
for k=1:K-1
    for i=1:N1*N2
        for j=1:N1*N2
            [rowl,coll]=ind2sub([N1 N2],i);
            [rowl2,coll2]=ind2sub([N1 N2],j);
            DRxSim(j,i,k)=sqrt((XRxSim(rowl)-XRxSim(rowl2))^2+(YRxSim(coll)-YRxSim(coll2))^2+(dzRx)^2);
            ThetaRxSim(j,i,k)=acos(dzRx/sqrt((XRxSim(rowl)-XRxSim(rowl2))^2+(YRxSim(coll)-YRxSim(coll2))^2+(dzRx)^2));
        end
    end
end

for k=1:K-1
    for i=1:N1*N2
        for j=1:N1*N2
            WRxSim(j,i,k)=dx*dy*cos(ThetaRxSim(j,i,k))/DRxSim(j,i,k)*(1/(2*pi*DRxSim(j,i,k))-1j/lambda)*exp(1j*2*pi*DRxSim(j,i,k)/lambda);
        end
    end
end


%%RxStream
for i=1:N1*N2
    for j=1:S
        [rowl,coll]=ind2sub([N1 N2],i);
        DRxStream(j,i)=sqrt((XRxSim(rowl)-XRxStream(j))^2+(YRxSim(coll)-YRxStream(j))^2+(ZRxSim(1)-ZRxStream(j))^2);
        ThetaRxStream(j,i)=acos(abs(ZRxSim(1)-ZRxStream(j))/sqrt((XRxSim(rowl)-XRxStream(j))^2+(YRxSim(coll)-YRxStream(j))^2+(ZRxSim(1)-ZRxStream(j))^2));
    end
end

for i=1:N1*N2
    for j=1:S
        WRxStream(j,i)=dx*dy*cos(ThetaRxStream(j,i))/DRxStream(j,i)*(1/(2*pi*DRxStream(j,i))-1j/lambda)*exp(1j*2*pi*DRxStream(j,i)/lambda);
    end
end




