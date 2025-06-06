clc
clear
close all


%%
for MonteCarloSeed=1:10 %Monte Carlo
    clearvars -except SEOptimize EEOptimize NMSEOptimize MonteCarloSeed SEOptimizeUB EEOptimizeUB
    disp(['**************************** Monte Carlo = ',num2str(MonteCarloSeed),' ****************************'])
    rng(MonteCarloSeed)
    Parameter
    TX_SIM_Initialize
    RX_SIM_Initialize
    Channel
    [E,Lam,F]=svd(DpG);
    lam=Lam(1:2*S,1:2*S);

    %%
    WaterFilling

    %%Initialize
    for seedInitialize=1:100
        rng(seedInitialize)
        Phi_Initialize
        TX_SIM
        RX_SIM
        DpH=DpGRxSim*DpG*DpGTxSim;
        alpha=inv(DpH(:)'*DpH(:))*DpH(:)'*lam(:);
        NMSEInitialize(seedInitialize)=norm(alpha*DpH-lam,'fro')^2/norm(lam,'fro')^2;
    end

    %%
    [~,index]=min(NMSEInitialize);
    rng(index)
    Phi_Initialize
    TX_SIM
    RX_SIM
    DpH=DpGRxSim*DpG*DpGTxSim;
    alpha=inv(DpH(:)'*DpH(:))*DpH(:)'*lam(:);
    InitializeDisp

    %%LGD-WF
    for round=2:Emax+1
        disp(['-----------------','Iteration round=',num2str(round-1),'-----------------'])
        UpdateTxPhi
        UpdateRxPhi
        Calculate
        eata=eata*beta;
    end

    SEOptimizeUB(MonteCarloSeed)=sum(log2(1+Ps.*diag(lam).^2/noise));
    EEOptimizeUB(MonteCarloSeed)=SEOptimizeUB(MonteCarloSeed)/P;
    SEOptimize(MonteCarloSeed)=SE(end);
    EEOptimize(MonteCarloSeed)=EE(end);
    NMSEOptimize(MonteCarloSeed)=NMSE(end);
    disp(['^^^^^^^^^^^^^^^^^^^^^^^^^^^ Optimized result MonteCarloSeed = ',num2str(MonteCarloSeed),' ^^^^^^^^^^^^^^^^^^^^^^^^^^^'])
    disp(['NMSEOptimize=',num2str(NMSEOptimize(MonteCarloSeed))])
    disp(['SEOptimize=',num2str(SEOptimize(MonteCarloSeed))])
    disp(['EEOptimize=',num2str(EEOptimize(MonteCarloSeed))])
    disp(['SEOptimizeUB=',num2str(SEOptimizeUB(MonteCarloSeed))])
    disp(['EEOptimizeUB=',num2str(EEOptimizeUB(MonteCarloSeed))])
    fprintf('\n\n\n\n\n')

end

%%
disp(['!!!!!!!!!!!!!!!!!!!!!!!!! Monte Carlo result !!!!!!!!!!!!!!!!!!!!!!!!!'])
SEOptimizeMonteCarlo=mean(SEOptimize);
EEOptimizeMonteCarlo=mean(EEOptimize);
SEOptimizeUBMonteCarlo=mean(SEOptimizeUB);
EEOptimizeUBMonteCarlo=mean(EEOptimizeUB);
NMSEOptimizeMonteCarlo=mean(NMSEOptimize);
disp(['NMSEOptimizeMonteCarlo=',num2str(NMSEOptimizeMonteCarlo)])
disp(['SEOptimizeMonteCarlo=',num2str(SEOptimizeMonteCarlo)])
disp(['EEOptimizeMonteCarlo=',num2str(EEOptimizeMonteCarlo)])
disp(['SEOptimizeUBMonteCarlo=',num2str(SEOptimizeUBMonteCarlo)])
disp(['EEOptimizeUBMonteCarlo=',num2str(EEOptimizeUBMonteCarlo)])

%%save
folderName = 'Result';
fileName = ['DPSIM-HMIMO,M=',num2str(M1*M2),',N=',num2str(N1*N2),',L=',num2str(L),',K=',num2str(K),',epsilon=',num2str(epsilon),',PdBm=',num2str(PdBm),'.mat'];
fullPath = fullfile(folderName, fileName);
save(fullPath,'NMSEOptimizeMonteCarlo','SEOptimizeMonteCarlo', 'EEOptimizeMonteCarlo','SEOptimizeUBMonteCarlo','EEOptimizeUBMonteCarlo');