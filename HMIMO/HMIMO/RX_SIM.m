
%%GRxSim
GRxSim=zeros(N1*N2,N1*N2);
GRxSim=PhiRxSim(:,:,1);
for k=1:K-1
    GRxSim=GRxSim*WRxSim(:,:,k)*PhiRxSim(:,:,k+1);
end
GRxSim=WRxStream*GRxSim;