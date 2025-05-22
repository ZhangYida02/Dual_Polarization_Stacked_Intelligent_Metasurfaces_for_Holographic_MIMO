%%
PhiRxSimPartDeri=zeros(N1*N2,K);
PhiRxSimPartDeriNorm=zeros(N1*N2,K);
for k=1:K
    for n=1:N1*N2
        for s1=1:S
            for s2=1:S

                if K==1
                    Y(s1,s2)=WRxStream(s1,n)*G(n,:)*GTxSim(:,s2);
                else
                if k==1
                    YLeft=WRxStream(s1,n);
                    YRight=G*GTxSim(:,s2);
                elseif k==K
                    YLeft=WRxStream(s1,:);
                    YRight=G(n,:)*GTxSim(:,s2);
                else
                    YLeft=WRxStream(s1,:);
                    YRight=G*GTxSim(:,s2);
                end

                %%
                for k1=1:1:k-1
                    if k1==k-1
                        YLeft=YLeft*PhiRxSim(:,:,k1)*WRxSim(:,n,k1);
                    else
                        YLeft=YLeft*PhiRxSim(:,:,k1)*WRxSim(:,:,k1);
                    end
                end
                for k2=K:-1:k+1
                    if k2==k+1
                        YRight=WRxSim(n,:,k2-1)*PhiRxSim(:,:,k2)*YRight;
                    else
                        YRight=WRxSim(:,:,k2-1)*PhiRxSim(:,:,k2)*YRight;
                    end
                end
                Y(s1,s2)=YLeft*YRight;
                end

            end
        end
        tempRx=imag(conj(alpha*PhiRxSim(n,n,k)*Y).*(alpha*H-lam));
        PhiRxSimPartDeri(n,k)=2*sum(sum(tempRx));
    end
    PhiRxSimPartDeriNorm(:,k)=PhiRxSimPartDeri(:,k)*pi/max(abs(PhiRxSimPartDeri(:,k)));
    PhiRxSim(:,:,k)=diag(exp(1j*(angle(diag(PhiRxSim(:,:,k)))-eata*PhiRxSimPartDeriNorm(:,k))));

    %%
    RX_SIM
    H=GRxSim*G*GTxSim;
    alpha=inv(H(:)'*H(:))*H(:)'*lam(:);
end

%%
