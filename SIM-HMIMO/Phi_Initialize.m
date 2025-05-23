%%PhiTxSimInit
PhiTxSimInit=zeros(M1*M2,M1*M2,L);
for l=1:L
    PhiTxSim(:,:,l)=diag(exp(1j*unifrnd(0,2*pi,[M1*M2 1])));
end

%%PhiRxSimInit
PhiRxSimInit=zeros(N1*N2,N1*N2,K);
for l=1:K
    PhiRxSim(:,:,l)=diag(exp(1j*unifrnd(0,2*pi,[N1*N2 1])));
end