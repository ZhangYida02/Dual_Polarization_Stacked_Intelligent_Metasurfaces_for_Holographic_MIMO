%%
DpPhiRxSimPartDeri=zeros(2*N1*N2,K);
DpPhiRxSimPartDeriNorm=zeros(2*N1*N2,K);
for k=1:K
    for n=1:2*N1*N2
        for s1=1:2*S
            for s2=1:2*S

                if L==1
                    Y(s1,s2)=DpWRxStream(s1,n)*DpG(n,:)*DpGTxSim(:,s2);
                else
                    if k==1
                        YLeft=DpWRxStream(s1,n);
                        YRight=DpG*DpGTxSim(:,s2);
                    elseif k==K
                        YLeft=DpWRxStream(s1,:);
                        YRight=DpG(n,:)*DpGTxSim(:,s2);
                    else
                        YLeft=DpWRxStream(s1,:);
                        YRight=DpG*DpGTxSim(:,s2);
                    end

                    %%
                    for k1=1:1:k-1
                        if k1==k-1
                            YLeft=YLeft*DpPhiRxSim(:,:,k1)*DpWRxSim(:,n,k1);
                        else
                            YLeft=YLeft*DpPhiRxSim(:,:,k1)*DpWRxSim(:,:,k1);
                        end
                    end
                    for k2=K:-1:k+1
                        if k2==k+1
                            YRight=DpWRxSim(n,:,k2-1)*DpPhiRxSim(:,:,k2)*YRight;
                        else
                            YRight=DpWRxSim(:,:,k2-1)*DpPhiRxSim(:,:,k2)*YRight;
                        end
                    end
                    Y(s1,s2)=YLeft*YRight;
                end

            end
        end
        tempRx=imag(conj(alpha*DpPhiRxSim(n,n,k)*Y).*(alpha*DpH-lam));
        DpPhiRxSimPartDeri(n,k)=2*sum(sum(tempRx));
    end
    DpPhiRxSimPartDeriNorm(:,k)=DpPhiRxSimPartDeri(:,k)*pi/max(abs(DpPhiRxSimPartDeri(:,k)));
    DpPhiRxSim(:,:,k)=diag(exp(1j*(angle(diag(DpPhiRxSim(:,:,k)))-eata*DpPhiRxSimPartDeriNorm(:,k))));

    %%
    RX_SIM
    DpH=DpGRxSim*DpG*DpGTxSim;
    alpha=inv(DpH(:)'*DpH(:))*DpH(:)'*lam(:);
end

%%
