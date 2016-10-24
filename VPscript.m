%function SD = VPmodelv1()
%%
r = 3e-9; % m, dDNA radius
Np = 3.5E16; %DNA/m2 x 10^4 from cm-2 unit

%w = 4E16; % electrostatic interaction parameter J.m2/mol for 1M Na+
As = 2*(3^0.5)*(r^2); %hexagon area
p = As*Np; %0.9
PM1 = getPM(p);
NM1 = getNM(As);

%%
%hyb result in 250k to 350k
SD = [];
QM = [];
for x = 280:2:380;
    QM1= getQM(NM1, x); 
    QMt = QM1';
    sd = PM1*QMt;
    SD = [SD sd];
    QM = [QM; QM1];
end
X1 = 280:2:380;
SDr= real(SD);
figure
plot(X1,SDr);
xlabel('Temp')
ylabel('Hyb ratio')
%%





