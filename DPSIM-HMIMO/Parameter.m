
%%EM
c=3*10^8;
f=28*10^9;
lambda=c/f;

%%TX-SIM
M1=10;
M2=10;
L=3;
DTx=0.05;
zTx=0;

%%RX-SIM
N1=10;
N2=10;
K=3;
DRx=0.05;
zRx=250;

%%stream
S=3;

%%noise
noisedB=-110;
noise=db2pow(noisedB)*10^-3;

%%P
PdBm=-10;
P=db2pow(PdBm)*10^-3;

%%algorithm
alpha=1;
eata=0.1;
beta=0.5;
epsilon=0.2;
Emax=20;

%%SIM
dx=lambda/2;
dy=lambda/2;
dzTx=DTx/L;
dzRx=DRx/K;


