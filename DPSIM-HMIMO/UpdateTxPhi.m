%%
DpPhiTxSimPartDeri=zeros(2*M1*M2,L);
DpPhiTxSimPartDeriNorm=zeros(2*M1*M2,L);
for l=1:L
    for m=1:2*M1*M2
        for s1=1:2*S
            for s2=1:2*S

                %%
                if L==1
                    X(s1,s2)=DpGRxSim(s1,:)*DpG(:,m)*DpWTxStream(m,s2);
                else
                    if l==1
                        XLeft=DpGRxSim(s1,:)*DpG;
                        XRight=DpWTxStream(m,s2);
                    elseif l==L
                        XLeft=DpGRxSim(s1,:)*DpG(:,m);
                        XRight=DpWTxStream(:,s2);
                    else
                        XLeft=DpGRxSim(s1,:)*DpG;
                        XRight=DpWTxStream(:,s2);
                    end

                    for l1=L:-1:l+1
                        if l1==l+1
                            XLeft=XLeft*DpPhiTxSim(:,:,l1)*DpWTxSim(:,m,l1-1);
                        else
                            XLeft=XLeft*DpPhiTxSim(:,:,l1)*DpWTxSim(:,:,l1-1);
                        end
                    end
                    for l2=1:1:l-1
                        if l2==l-1
                            XRight=DpWTxSim(m,:,l2)*DpPhiTxSim(:,:,l2)*XRight;
                        else
                            XRight=DpWTxSim(:,:,l2)*DpPhiTxSim(:,:,l2)*XRight;
                        end
                    end
                    X(s1,s2)=XLeft*XRight;
                end

            end
        end
        tempTx=imag(conj(alpha*DpPhiTxSim(m,m,l)*X).*(alpha*DpH-lam));
        DpPhiTxSimPartDeri(m,l)=2*sum(sum(tempTx));
    end
    DpPhiTxSimPartDeriNorm(:,l)=DpPhiTxSimPartDeri(:,l)*pi/max(abs(DpPhiTxSimPartDeri(:,l)));
    DpPhiTxSim(:,:,l)=diag(exp(1j*(angle(diag(DpPhiTxSim(:,:,l)))-eata*DpPhiTxSimPartDeriNorm(:,l))));

    %%
    TX_SIM
    DpH=DpGRxSim*DpG*DpGTxSim;
    alpha=inv(DpH(:)'*DpH(:))*DpH(:)'*lam(:);
    

end
