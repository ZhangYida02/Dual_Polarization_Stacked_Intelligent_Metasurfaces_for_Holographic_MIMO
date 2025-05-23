
%%DpGTxSim
DpGTxSim=zeros(2*M1*M2,2*M1*M2);

DpGTxSim=DpPhiTxSim(:,:,1);
for l=1:L-1
    DpGTxSim=DpPhiTxSim(:,:,l+1)*DpWTxSim(:,:,l)*DpGTxSim;
end
DpGTxSim=DpGTxSim*DpWTxStream;