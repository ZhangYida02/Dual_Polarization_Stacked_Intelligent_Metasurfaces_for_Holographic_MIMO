

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

SE(round)=capWF;
EE(round)=SE(round)/P;
NMSE(round)=norm(alpha*H-lam,'fro')^2/norm(lam,'fro')^2;

disp(['NMSE=',num2str(NMSE(round))])
disp(['SE=',num2str(SE(round))])
disp(['EE=',num2str(EE(round))])

