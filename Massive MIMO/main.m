clc
clear
close all


c=3*10^8;
f=28*10^9;
lambda=c/f;
zTx=0;
zRx=250;
noisedBm=-110;
noise=db2pow(noisedBm)*10^-3;
alp=0.2;
rx=32;
tx=256;
PL=20*log10(4*pi/lambda)+10*3.5*log10(zRx)+9;
rho_sq=1/db2pow(PL);

for PdBm=[-10:5:40]
    P=db2pow(PdBm)*10^-3;
    for seed =1:100
        rng(seed)
        G = sqrt((1-alp)*rho_sq/2) * (randn(rx,tx) + 1j*randn(rx,tx)); 
        SE(seed) = calculate_spectral_efficiency(G, P, noise);
        EE(seed) = SE(seed)/P;
    end
    disp(['PdBm=',num2str(PdBm),',massive MIMO SE=',num2str(mean(SE))])
    disp(['PdBm=',num2str(PdBm),',massive MIMO EE=',num2str(mean(EE))])
    SEMean(find(PdBm==[-10:5:40]))=mean(SE);
    EEMean(find(PdBm==[-10:5:40]))=mean(EE);
    disp('-----------------------------------------------------------')
end