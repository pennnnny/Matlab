%%
% Create a variable containing the location of the Xba array CEL files. For
% example:
Xba_celPath = 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\Hind\GSM248457.CEL';
%Xba_celPath = 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\Xba\100K_trios\chip_data\Xba';

%%
% Create a variable containing the location of the CDF library file for the
% EA 50KXba SNP array. For example:
libPath = 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\Hind\50k_hind240_libraryfile_rev3\CD_Mapping50K_Hind240_rev3\Full\Mapping50K_Hind240\LibFiles';
%xba 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\CD_Mapping50K_Xba240_rev3\Full\Mapping50K_Xba240\LibFiles';

%%
% |SCLC_Sample_CEL.txt|, a file provided with the Bioinformatics
% Toolbox(TM) software, contains a list of the 24 CEL file names used for
% this example, and the samples (5 SCLC cell lines and 19 primary tumors)
% to which they belong. Load this data into two MATLAB(R) variables. 
%fid = fopen('SCLC_Sample_CEL.txt','r');
%ftext = textscan(fid, '%q%q');
%fclose(fid);
%samples = ftext{1};
%cels = ftext{2};
%cels = 'Mapping50K_Xba240.CDF';
%nSample = numel(samples)

%% Accessing SNP Array Probe-Level Data
% The Affymetrix 50KXba SNP array has a density up to 50K SNP sites. Each
% SNP on the array is represented by a collection of probe quartets. A
% probe quartet consists of a set of probe pairs for both alleles (A and B)
% and for both forward and reverse strands (antisense and sense) for the
% SNP. Each probe pair consists a perfect match (PM) probe and a mismatch
% (MM) probe. The Bioinformatics Toolbox software provides functions to
% access the probe-level data.
%%
% The function |affyread| reads the CEL files and the CDF library files
% for Affymetrix SNP arrays. 
%%
% Read the sixth CEL file of the EA 50KXba data into a MATLAB structure.
%s_cel = affyread(fullfile(Xba_celPath, [cels{6} '.CEL']))
s_cel = affyread(fullfile(Xba_celPath, 'GSM248457.CEL')) %AD362T.CEL used in example
%%
% Read the CDF library file for the EA 50KXba array into a MATLAB
% structure.
s_cdf = affyread(fullfile(libPath, 'Mapping50K_Hind240.cdf'))
%%
% You can inspect the overall quality of the array by viewing the
% probe-level intensity data using the function |maimage|.
 maimage(s_cel)
%%
% The |affysnpquartets| function creates a table of probe quartets for a
% SNP. On Affymetrix 100K SNP arrays, a probe quartet contains 20 probe
% pairs. For example, to get detailed information on probe set number 6540,
% you can type the following commands:
ps_id = 6540; %6540/59024
ps_qt = affysnpquartets(s_cel, s_cdf, ps_id)

%%
% You can also view the heat map of the intensities of the PM and MM probe
% pairs of a SNP probe quartet using the |probesetplot| function. Click the
% *Insert Colorbar* button to show the color scale of the heat map.
probesetplot(s_cel, s_cdf, ps_id, 'imageonly', true);

%%
% In this view, the 20 probe pairs are ordered from left to right. The
% first two rows (10 probe pairs) correspond to allele A, and the last two
% rows (10 probe pairs) corresponds to allele B. For each allele, the left
% 5 probe pairs correspond to the sense strand (-), while the right 5 probe
% pairs correspond to the antisense (+) strand.

%% Importing and Converting the Data Set
% You will use the |celintensityread| function to read all 24 CEL files.
% The |celintensityread| function returns a structure containing the
% matrices of PM and MM (optional) intensities for the probes and their
% group numbers. In each probe intensity matrix, the column indices
% correspond to the order in which the CEL files were read, and each row
% corresponds to a probe. For copy number (CN) analysis, only PM
% intensities are needed.

%%
% Import the probe intensity data of all EA 50KXba arrays into a MATLAB
% structure.
FAL = logical(0);
XbaDataOg = celintensityread('GSM248457.CEL', 'Mapping50K_Hind240.cdf',...
                     'celpath', Xba_celPath, 'cdfpath', libPath, 'PMONLY', FAL)
                 
%% Probe-Level Summarization
% By using the |GroupNumbers| field data from the structure |XbaData|, you
% can extract the intensities for allele A and allele B for each probe. Use
% the function |affysnpintensitysplit| to split the probe intensities
% matrix |PMIntensities| into two single-precision matrices,
% |PMAIntensities| and |PMBIntensities|, for allele A and allele B probes
% respectively. The number of probes in each matrix is the maximum number
% of probes for each allele. 
XbaData = affysnpintensitysplit(XbaDataOg)

%%
% For total copy number analysis, a simplification is to ignore the allele
% A and allele B sequences and their strand information and, instead,
% combine the PM intensities for allele A and allele B of each probe pair. 
PM_Xba = XbaData.PMAIntensities;
%PM_XbaB = XbaData.PMBIntensities;% + XbaData.PMBIntensities;

%%
% For a particular SNP, we now have _K_ (_K=5_ for Affymetrix Mapping 100K
% arrays) added signals, each signal being a measure of the same thing
% - the total CN. However, each of the _K_ signals has slightly different
% sequences, so their hybridization efficiency might differ. You can use
% RMA summarization methods to sum up allele probe intensities for each SNP
% probe set.
PM_Xba = rmasummary(XbaData.ProbeIndices, XbaData.PMAIntensities);   

%%
%read probe seq
%S = affyprobeseqread('Mapping50K_Xba240_probe_fasta', 'Mapping50K_Xba240.CDF',...
%               'seqpath', seqPath, 'cdfpath', libPath, 'seqonly', true)
%% load probe sequence
% The function |affysnpannotread| reads the annotation file and returns a
% structure containing SNP chromosome information, chromosomal positions,
% sequences and PCR fragment length information ordered by probe set IDs
% from the second input variable.
annoPath = 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\Hind\Mapping50K_Hind240.na32.annot.csv';
%annoPath = 'C:\Users\HowardNYU\Documents\MATLAB\affymetrix_matlab_demo\Mapping100K.na32.annot.db';
annoFile = fullfile(annoPath, 'Mapping50K_Hind240.na32.annot.csv'); %xba 'Mapping50K_Xba240.na31.annot.csv'
annot_Xba = affysnpannotread(annoFile, XbaDataOg.ProbeSetIDs)

%%
%find match ProbesetID
Pind =zeros(length(XbaData.ProbeSetIDs),1);
for x = 1:length(XbaData.ProbeSetIDs)
    pind = find(strncmp(annot_Xba.ProbeSetIDs,XbaData.ProbeSetIDs{x}, length(XbaData.ProbeSetIDs{x})));
    if isempty(pind)
        pind = 1; %not found probesetID to make value as 0*** (should use other)
        Pind(x,1) = pind;
    else        
        Pind(x,1) = pind;
    end
end

%% remove mismatch 
%focus on allel A
sequence_Xba = annot_Xba.Sequence;
ss = size(sequence_Xba);
Sequence_Xba = cell(zeros);
Sequence_Xbaf = cell(zeros);
Sequence_Xbafr = cell(zeros);
Sequence_Xbab = cell(zeros);
Sequence_Xbabr = cell(zeros);
%Generate new array containing probe sequences
for i = 1: ss(1)
    sxba = sequence_Xba{i};
    Sequence_Xba = [Sequence_Xba; sxba];
end
% select allel +subfunction apartseq()
Countcheck = [];
for i  =1:ss(1)
    checksxb = find(Sequence_Xba{i} == '[');
    checknullsxb = isempty(checksxb);
    if checknullsxb 
        App = 'N'; %for allel A
        Appr = 'N';
        Bpp = 'N';
        Bppr = 'N';
        countcheck = 1;
    else
        [App Bpp Appr Bppr] = getProbeseq(Sequence_Xba{i}); %get repeat part of probe sequence
        countcheck = 0;
    end
    Sequence_Xbaf = [Sequence_Xbaf; App]; 
    Sequence_Xbafr = [Sequence_Xbafr; Appr]; 
    Sequence_Xbab = [Sequence_Xbab; Bpp]; 
    Sequence_Xbabr = [Sequence_Xbabr; Bppr]; 
    Countcheck = [Countcheck countcheck];
end
%%
% Run UNAfold to get Thermodynamic value
n = 200; %each batch max strand for UNAfold *****n>400 can't run**** %only 200 for hyb
%*** seq selector
SEQ = Sequence_Xbaf;
SEQr = Sequence_Xbafr;
%get number of submisson 
endloop = ceil((length(SEQ)-55)/n); %first 55 sequnce is blank****
tic
DG = cell(length(SEQ)-55, 1);
DH = cell(length(SEQ)-55, 1);
DS = cell(length(SEQ)-55, 1);
TM = cell(length(SEQ)-55, 1);
for x = 1:endloop
    if 55+n*x <= length(SEQ)
     	seq = [SEQ{56+n*(x-1):55+n*x}];
        seqr = [SEQr{56+n*(x-1):55+n*x}];
    else
        seq = [SEQ{56+n*(x-1):length(SEQ)}];
        seqr = [SEQr{56+n*(x-1):length(SEQr)}];
    end
    %seq < n 
    if 55+n*x <= length(SEQ)
        [dg dh ds tm] = getWebthermo(seq, seqr);
        DG(1+n*(x-1):n*x, 1) = dg;
        DH(1+n*(x-1):n*x, 1) = dh;
        DS(1+n*(x-1):n*x, 1) = ds;
        TM(1+n*(x-1):n*x, 1) = tm;
    else
        [dg dh ds tm] = getWebthermo(seq, seqr);
        DG(1+n*(x-1):length(DG), 1) = dg;
        DH(1+n*(x-1):length(DH), 1) = dh;
        DS(1+n*(x-1):length(DS), 1) = ds;
        TM(1+n*(x-1):length(TM), 1) = tm;
    end
    
end 
toc
%%
%********************UNAfold %converter 
%get thermo value from UNAfold cell
DGbhybd = zeros(length(DGbhyb),1);
for x =1:length(DGbhyb)
    dg37 = DGbhyb{x};
    dg37d = str2double(dg37);
    DGbhybd(x,1) = dg37d;
end