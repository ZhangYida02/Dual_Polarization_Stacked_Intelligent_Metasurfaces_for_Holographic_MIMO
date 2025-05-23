
%%GRxSim
DpGRxSim=zeros(2*N1*N2,2*N1*N2);
DpGRxSim=DpPhiRxSim(:,:,1);
for k=1:K-1
    DpGRxSim=DpGRxSim*DpWRxSim(:,:,k)*DpPhiRxSim(:,:,k+1);
end
DpGRxSim=DpWRxStream*DpGRxSim;