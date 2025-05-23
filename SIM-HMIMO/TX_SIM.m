
%%GTxSim
GTxSim=zeros(M1*M2,M1*M2);

GTxSim=PhiTxSim(:,:,1);
for l=1:L-1
    GTxSim=PhiTxSim(:,:,l+1)*WTxSim(:,:,l)*GTxSim;
end
GTxSim=GTxSim*WTxStream;