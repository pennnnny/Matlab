%%
%get index for XbaData.ProbeIndices
XbaSeninfo = zeros(length(XbaData.ProbeSetIDs),10);
for x = 1:length(XbaData.ProbeSetIDs)
    pindID = [1+((x-1)*10):(1+((x-1)*10)+9)];
    XbaSeninfo(x,:) = pindID;
end
mIntA = zeros(length(s_cdf.ProbeSets),1); % = PMAmean
mIntAs = zeros(length(s_cdf.ProbeSets),1); % = iPMAsen1
mIntAas = zeros(length(s_cdf.ProbeSets),1); % = iPMAsen2

mIntB = zeros(length(s_cdf.ProbeSets),1); % = PMBmean
mIntBs = zeros(length(s_cdf.ProbeSets),1); % = iPMBsen1
mIntBas = zeros(length(s_cdf.ProbeSets),1); % = iPMBsen2
pA = [1,3,5,7,9,11,13,15,17,19]; %Allele A probe index
pB = [2,4,6,8,10,12,14,16,18,20];   %Allele B probe index
for i = 1:length(s_cdf.ProbeSets)
    pn = s_cdf.ProbeSets(i).Name;
    if ~isempty(pn)
        iXba = find(strncmp(pn, XbaData.ProbeSetIDs, length(pn)));
        if ~isempty(iXba)
            seninfo = s_cdf.ProbeSets(i).ProbePairs(:,2); %sense id
            seninfoA = seninfo(pA);
            seninfoB = seninfo(pB);
            %read A and B intensity
            intIndx = XbaSeninfo(iXba,:);
            intA = XbaData.PMAIntensities(intIndx); %intensity A
            intB = XbaData.PMBIntensities(intIndx);
            mintA = mean(intA);
            mintB = mean(intB);
            %read sense and antisense seperately
            senix = find(seninfoA == 1); 
            intAs = intA(senix);
            mintAs = mean(intAs);
            antisenix = find(seninfoA == 2); 
            intAas = intA(antisenix);
            mintAas = mean(intAas);
            
            senixb = find(seninfoB == 1); 
            intBs = intB(senixb);
            mintBs = mean(intBs);
            antisenixb = find(seninfoB == 2); 
            intBas = intB(antisenixb);
            mintBas = mean(intBas);
        else
            mintA = 0;
            mintAs = 0;
            mintAas = 0;
            mintB = 0;
            mintBs = 0;
            mintBas = 0;
        end
    else
        mintA = 0;
        mintAs = 0;
        mintAas = 0;
        mintB = 0;
        mintBs = 0;
        mintBas = 0;
    end    
    mIntA(i) = mintA;
    mIntAs(i) = mintAs;
    mIntAas(i) = mintAas;
    mIntB(i) = mintB;
    mIntBs(i) = mintBs;
    mIntBas(i) = mintBas;
end
PMAmean = mIntA(56:57299); %get non zero values %Xba 59015 / Hind = 57299
PMBmean = mIntB(56:57299);
iPMAsen1 = mIntAs(56:57299);
iPMAsen2 = mIntAas(56:57299);
iPMBsen1 = mIntBs(56:57299);
iPMBsen2 = mIntBas(56:57299);

%remove B is stroger
RPM = PMBmean./PMAmean; 
RPMr = PMAmean./PMBmean; 
tiPM = zeros(length(iPMAsen1),8);
%tiPMAsen2 = zeros(length(iPMAsen2),1);
for x = 1:length(iPMAsen1)
    if RPM(x) <= 0.5823  %SNP in Allel AA***** % RPM(x) <= 0.5816 //// DPM(x) >= 0 
        %if std(pmasumsen1(x,:)) < mean(pmasumsen1(x,:))*0.05 
        tpma = iPMAsen1(x);
        tpma2 = iPMAsen2(x);
        tpmb = 0;
        tpmb2 = 0;
        tpmaba = 0;
        tpmaba2 = 0;
        tpmabb = 0;
        tpmabb2 = 0;
    elseif RPMr(x) <= 0.5823 %SNP = BB
        tpma = 0;
        tpma2 = 0;
        tpmb = iPMBsen1(x);
        tpmb2 = iPMBsen2(x);
        tpmaba = 0;
        tpmaba2 = 0;
        tpmabb = 0;
        tpmabb2 = 0;
    elseif RPM(x) > 0.5823 && RPMr(x) > 0.5823 % SNP =AB
        tpma = 0;
        tpma2 = 0;
        tpmb = 0;
        tpmb2 = 0;
        tpmaba = iPMAsen1(x);
        tpmaba2 = iPMAsen2(x);
        tpmabb = iPMBsen1(x);
        tpmabb2 = iPMBsen2(x);
    else
        tpma = 0;
        tpma2 = 0;
        tpmb = 0;
        tpmb2 = 0;
        tpmaba = 0;
        tpmaba2 = 0;
        tpmabb = 0;
        tpmabb2 = 0;   
            
    end
    tiPM(x,1:8) = [tpma, tpma2, tpmb, tpmb2, tpmaba, tpmaba2, tpmabb, tpmabb2];
    %mean signal for each probe sorted by probe sequences
    %(A/B/sense/antisense)
end
%get mean/median at each dG with N data points
% [A A B B A A B B] 
%DG45hybd8 = [DG45hybd,DG45hybd,DGbhybd,DGbhybd,DG45hybd,DG45hybd,DGbhybd,DGbhybd];
DG45hybd8 = [DGAhybd,DGAhybd,DGBhybd,DGBhybd,DGAhybd,DGAhybd,DGBhybd,DGBhybd];
dgh = [-5:-0.1:-20];
tiPMg = zeros(length(dgh),1);
tiSen2h = zeros(length(dgh),1);
Ndot = zeros(length(dgh),1);
STD = zeros(length(dgh),1);
Nbar = 100;
for x = 1:length(dgh)
    g = dgh(x); %find value for each dG
    indx = find(DG45hybd8 == g);
    tipmg = tiPM(indx);
    nz = find(tipmg ~= 0);
    Indxnz = indx(nz);
    ndot = length(Indxnz);
        if length(Indxnz) >= Nbar
            isen1 = mean(tiPM(Indxnz)); %%**** check value
            st = std(tiPM(Indxnz));
        else
            isen1 = 0;
            st = 0;
        end    
    tiPMg(x) = isen1;
    Ndot(x) = ndot;
    STD(x) = st;
end
figure
plot(dgh, tiPMg, '.')
%%
%histogram plot of AA AB BB
rPM = zeros(length(PMAmean),1);
for x = 1:length(PMAmean)
    if PMAmean(x)-PMBmean(x) >=0 % allel AA or AB
        rpm = PMBmean(x)/PMAmean(x);
    else
        rpm = PMAmean(x)/PMBmean(x);
    end
    rPM(x) = rpm;
end
%figure;plot([1:58991], rPM, '.')
[fr,xr] = ksdensity(rPM);
figure;plot(xr, fr/sum(fr))
nr0 = histc(rPM, [0:0.01:1]);
figure;bar([0:0.01:1],nr0,'histc')
%***** detemine cutoff value of SNP type and rerun above section******
%%
% check the effect of pararmeters 
%get signal with sequence info
pmasort = sort(PMAmean);
ib10 = pmasort(length(PMAmean)-1000:length(PMAmean));
id10 = pmasort(1:1000);
mn =[1,2,4,6,8,10,20,30,40,50,60,70,80,90,100,200];
figure
hold on
%for m = 1:length(mn)
I02 = (mean(ib10) - mean(id10));  %****mn(m); %second parameter
Ibg2 = mean(id10); %first parameter
c12 = 45; % target = 100uM (thrid parameter)
T = 318.15; % for probe

DDG = zeros(length(iPMAsen1),8);
Term1 = zeros(length(iPMAsen1),8);
Term2 = zeros(length(iPMAsen1),8);
Term3 = zeros(length(iPMAsen1),8);
    for i = 1:8
        tiPMi = tiPM(:,i);
        for x = 1:length(tiPMi)
            if tiPMi(x) ~=0
                ddgc = ((8.314)*T*log(log(I02/(tiPMi(x)-Ibg2)))-(8.314*T)*log(c12))-DG45hybd8(x,i)*1000; %-DG45hybd8(x,i)*1000
                term1 = (8.314)*T*log(log(I02/(tiPMi(x)-Ibg2)));
                term2 = (8.314*T)*log(c12);
                term3 = DG45hybd8(x,i)*1000;
                if isreal(ddgc)
                    ddg = ddgc;
                else
                    ddg = 0;
                end
            else
                ddg = 0;
                term1 = 0;
                term2 = 0;
                term3 = 0;
            end
            DDG(x,i) = ddg;
            Term1(x,i) = term1;
            Term2(x,i) = term2;
            Term3(x,i) = term3;
        end
    end
    %got DDG for each sequence

    % plot dgpp for each value
    dg  =[3:-0.1:-10]; %[3:-0.1:0] /[-0.1:-0.1:-10]
    dpp = round(dg.*10)./10;
    %[As Aas Bs Bas As Aas Bs Bas]
    DG45sd8 = [DGAsd,DGAasd,DGBsd,DGBasd,DGAsd,DGAasd,DGBsd,DGBasd];
    DG8 = round(DG45sd8.*10)./10;
    mDDG = zeros(length(dg),1);
    mTerm1 = zeros(length(dg),1);
    NdotDDG = zeros(length(dg),1);
    StdDDG = zeros(length(dg),1);
    ER = zeros(length(dg),2);
    ERt = zeros(length(dg),2);
    Nbar = 100;
    for x = 1:length(dg)
        dgx = dpp(x); %find value for each dG
        indx = find(DG8 == dgx);
        nz = find(DDG(indx) ~= 0);
        indxnz = indx(nz);
        ndotg = length(nz);
        stdh = std(DDG(indxnz));
        DDGi = DDG(indxnz);
        Term1i = Term3(indxnz); % ***** change Term 
        %DGsurf = Term1(indxnz)-Term2(indxnz);/DGsol = Term3(indxnz)
            if length(indxnz) >= Nbar
                isen1 = mean(DDGi); %%**** check value
                term1i = mean(Term1i);
            else
                isen1 = 0;
                term1i = 0;
            end
            %get square root erro for both upper and lower value
                Eup = zeros(length(DDGi),1);
                Elow = zeros(length(DDGi),1);
                counteup = 0;
                countelow = 0;
            if length(indxnz) >= Nbar
                for i = 1:length(indxnz)
                    if isen1 > DDGi(i) %get upper error
                        deup = (DDGi(i) - isen1).^2;
                        delow = 0;
                        counteup = 1 + counteup;
                    else
                        delow = (isen1 - DDGi(i)).^2;
                        deup = 0;
                        countelow = 1 + countelow;
                    end
                    Eup(i) = deup;
                    Elow(i) = delow;
                end
                mEup = (sum(Eup)/counteup)^0.5;
                mElow = (sum(Elow)/countelow)^0.5;
            else
                mEup = 0;
                mElow = 0;
            end
            % get error for term
            if length(indxnz) >= Nbar
                for i = 1:length(indxnz)
                    if term1i > Term1i(i) %get upper error
                        deup = (Term1i(i) - term1i).^2;
                        delow = 0;
                        counteup = 1 + counteup;
                    else
                        delow = (term1i - Term1i(i)).^2;
                        deup = 0;
                        countelow = 1 + countelow;
                    end
                    Eup(i) = deup;
                    Elow(i) = delow;
                end
                mEupt = (sum(Eup)/counteup)^0.5;
                mElowt = (sum(Elow)/countelow)^0.5;
            else
                mEupt = 0;
                mElowt = 0;
            end

        ER(x,1) = mEup;
        ER(x,2) = mElow;
        ERt(x,1) = mEupt;
        ERt(x,2) = mElowt;
        mDDG(x) = isen1;
        mTerm1(x) = term1i;
        NdotDDG(x) = ndotg;
        StdDDG(x) = stdh;
    end
    inz = find(mDDG ~= 0);
    plot(dg(inz), mDDG(inz), '.')
    %plot(dg, mTerm1, '.')
    
%end
hold off
%%
%get DGpp profile
ixnz = find(tiPM ~= 0);
ppnz = DG45sd8(ixnz);
nr = histc(ppnz, [-10:0.01:3]);
figure;bar( [-10:0.01:3],nr,'histc')
