NMSEInitialize=norm(alpha*H-lam,'fro')^2/norm(lam,'fro')^2;
disp(['-----------------Initialize-----------------'])
disp(['NMSE=',num2str(NMSEInitialize)])
NMSE(1)=NMSEInitialize;
capWF=0;
for s1=1:S
    inteWF=0;
    for s2=1:S
        if s2~=s1
            inteWF=inteWF+Ps(s2)*abs(alpha*H(s1,s2)).^2;
        end
    end
    capWF=capWF+log2(1+Ps(s1)*abs(alpha*H(s1,s1)).^2/(inteWF+noise));
end
SE(1)=capWF;
EE(1)=SE(1)/P;
disp(['SE=',num2str(SE(1))])
disp(['EE=',num2str(EE(1))])