%%PhiTxSimInit
DpPhiTxSim=zeros(2*M1*M2,2*M1*M2,L);
for l=1:L
    DpPhiTxSim(:,:,l)=diag(exp(1j*unifrnd(0,2*pi,[2*M1*M2 1])));
end

%%PhiRxSimInit
DpPhiRxSim=zeros(2*N1*N2,2*N1*N2,K);
for l=1:K
    DpPhiRxSim(:,:,l)=diag(exp(1j*unifrnd(0,2*pi,[2*N1*N2 1])));
end