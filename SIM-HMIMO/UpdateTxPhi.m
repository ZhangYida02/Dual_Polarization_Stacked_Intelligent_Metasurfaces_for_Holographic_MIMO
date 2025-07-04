%%
PhiTxSimPartDeri=zeros(M1*M2,L);
PhiTxSimPartDeriNorm=zeros(M1*M2,L);
for l=1:L
    for m=1:M1*M2
        for s1=1:S
            for s2=1:S

                %%
                if L==1
                    X(s1,s2)=GRxSim(s1,:)*G(:,m)*WTxStream(m,s2);
                else

                    if l==1
                        XLeft=GRxSim(s1,:)*G;
                        XRight=WTxStream(m,s2);
                    elseif l==L
                        XLeft=GRxSim(s1,:)*G(:,m);
                        XRight=WTxStream(:,s2);
                    else
                        XLeft=GRxSim(s1,:)*G;
                        XRight=WTxStream(:,s2);
                    end

                    for l1=L:-1:l+1
                        if l1==l+1
                            XLeft=XLeft*PhiTxSim(:,:,l1)*WTxSim(:,m,l1-1);
                        else
                            XLeft=XLeft*PhiTxSim(:,:,l1)*WTxSim(:,:,l1-1);
                        end
                    end
                    for l2=1:1:l-1
                        if l2==l-1
                            XRight=WTxSim(m,:,l2)*PhiTxSim(:,:,l2)*XRight;
                        else
                            XRight=WTxSim(:,:,l2)*PhiTxSim(:,:,l2)*XRight;
                        end
                    end
                    X(s1,s2)=XLeft*XRight;
                end

            end
        end
        tempTx=imag(conj(alpha*PhiTxSim(m,m,l)*X).*(alpha*H-lam));
        PhiTxSimPartDeri(m,l)=2*sum(sum(tempTx));
    end
    PhiTxSimPartDeriNorm(:,l)=PhiTxSimPartDeri(:,l)*pi/max(abs(PhiTxSimPartDeri(:,l)));
    PhiTxSim(:,:,l)=diag(exp(1j*(angle(diag(PhiTxSim(:,:,l)))-eata*PhiTxSimPartDeriNorm(:,l))));

    %%
    TX_SIM
    H=GRxSim*G*GTxSim;
    alpha=inv(H(:)'*H(:))*H(:)'*lam(:);   

end
