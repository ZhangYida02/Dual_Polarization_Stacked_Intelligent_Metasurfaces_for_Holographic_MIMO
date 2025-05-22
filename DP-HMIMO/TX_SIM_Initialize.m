
%%TX_SIM
XTxSim=-lambda/2*(M1-1)/2:lambda/2:lambda/2*(M1-1)/2;
YTxSim=-lambda/2*(M2-1)/2:lambda/2:lambda/2*(M2-1)/2;
ZTxSim=0:dzTx:dzTx*(L-1);
WTxSim=zeros(M1*M2,M1*M2,L-1);
DpWTxSim=zeros(2*M1*M2,2*M1*M2,L-1);
DTxSim=zeros(M1*M2,M1*M2,L-1);
ThetaTxSim=zeros(M1*M2,M1*M2,L-1);


%%TxStream
DTxStream=zeros(M1*M2,S);
WTxStream=zeros(M1*M2,S);
DpWTxStream=zeros(2*M1*M2,2*S);
ThetaTxStream=zeros(M1*M2,S);
XTxStream=-lambda/2*(S-1)/2:lambda/2:lambda/2*(S-1)/2;
YTxStream=zeros(1,S);
ZTxStream=-dzTx*ones(1,S);
[X_gird,Y_gird]=meshgrid(XTxStream,YTxStream);


%%TX_SIM
for l=1:L-1
    for i=1:M1*M2
        for j=1:M1*M2
            [rowl,coll]=ind2sub([M1 M2],i);
            [rowl2,coll2]=ind2sub([M1 M2],j);
            DTxSim(j,i,l)=sqrt((XTxSim(rowl)-XTxSim(rowl2))^2+(YTxSim(coll)-YTxSim(coll2))^2+(dzTx)^2);
            ThetaTxSim(j,i,l)=acos(dzTx/sqrt((XTxSim(rowl)-XTxSim(rowl2))^2+(YTxSim(coll)-YTxSim(coll2))^2+(dzTx)^2));
        end
    end
end

for l=1:L-1
    for i=1:M1*M2
        for j=1:M1*M2
            WTxSim(j,i,l)=dx*dy*cos(ThetaTxSim(j,i,l))/DTxSim(j,i,l)*(1/(2*pi*DTxSim(j,i,l))-1j/lambda)*exp(1j*2*pi*DTxSim(j,i,l)/lambda);
        end
    end
end


%%TxStream
for i=1:S
    for j=1:M1*M2
        [rowl,coll]=ind2sub([M1 M2],j);
        DTxStream(j,i)=sqrt((XTxSim(rowl)-XTxStream(i))^2+(YTxSim(coll)-YTxStream(i))^2+(ZTxSim(1)-ZTxStream(i))^2);
        ThetaTxStream(j,i)=acos(abs(ZTxSim(1)-ZTxStream(i))/sqrt((XTxSim(rowl)-XTxStream(i))^2+(YTxSim(coll)-YTxStream(i))^2+(ZTxSim(1)-ZTxStream(i))^2));
    end
end

for i=1:S
    for j=1:M1*M2
        WTxStream(j,i)=dx*dy*cos(ThetaTxStream(j,i))/DTxStream(j,i)*(1/(2*pi*DTxStream(j,i))-1j/lambda)*exp(1j*2*pi*DTxStream(j,i)/lambda);
    end
end


%%DP
for l=1:L-1
    DpWTxSim(:,:,l)=[WTxSim(:,:,l) zeros(size(WTxSim(:,:,l)));zeros(size(WTxSim(:,:,l))) WTxSim(:,:,l)];
end

DpWTxStream=[WTxStream zeros(size(WTxStream));zeros(size(WTxStream)) WTxStream];


