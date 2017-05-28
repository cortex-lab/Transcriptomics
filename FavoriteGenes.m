function GeneList = FavoriteGenes();
GeneList = {'Adarb2', 'Adcy1', 'Akr1c18', ...
    'Cacna2d1', 'Cacna2d2', 'Cadps2', 'Calb1', 'Calb2', 'Caln1', 'Car2', 'Cbln2', 'Cck', 'Cd24a', 'Cdh13', 'Cend1', ...
    'Chrm2', 'Cidea', ... 
    'Cnr1', 'Col19a1', 'Cox6a2', 'Cplx2', 'Cplx3', 'Cpne5', 'Crabp1', 'Crh', 'Crhbp', 'Cryab', 'Csrp2', ...
    'Cxcl14', 'Dpy19l1', 'Elfn1', 'Enpp2', 'Fxyd6', 'Gabra1', 'Gabrd', ...
    'Gria3', 'Grik1', 'Grin3a', 'Grm1', ...
    'Hapln1', 'Hpcal1', 'Htr2c', 'Htr3a', 'Id2', 'Igf1', 'Igfbp5', 'Kcnc1', 'Kcnip2', 'Kctd12', 'Kit', 'Krt73', ...
    'Lamp5', 'Lgals1', 'Lhx6', ...
    'Lypd6','Lypd6b', 'Mafb', 'Mef2c', 'Ndnf', 'Necab1', 'Necab2', ...
    'Nfib', 'Nos1', 'Npas1', 'Npas3', 'Npy', 'Npy1r', 'Npy2r', 'Nr2f2',  'Nrsn1', 'Nrxn3', 'Nsg1', 'Ntng1', 'Nxph1', ...
    'Olfm3', 'Pcp4', 'Penk', 'Plcxd3', 'Pnoc', 'Ppargc1a', 'Prkacb', 'Pthlh', 'Ptprd', 'Ptpre', 'Pvalb', 'Qrfpr' ...
    'Rab3b', 'Rbp4', 'Reln', 'Rgs10', 'Rgs12', 'Rgs17', 'Rgs8', ...
    'Satb1', 'Sel1l3', 'Sema3c', 'Sema5a', 'Serpine2',  ...
    'Slc17a8', 'Slc7a14', 'Snca', 'Sncg', 'Stmn2', ...
    'Sox5', 'Sox6', 'Sparcl1', 'Spock1', 'Sst', 'Sstr1', 'Stk32a', 'Synpr', 'Syt2', 'Syt6', 'Tac1', 'Tac2', 'Tacr1', 'Thrsp', 'Thsd7a'...
    'Tnfaip8l3', 'Trp53i11', 'Trpc6', 'Tubb2a', 'Vip', 'Vwa5a', 'Wnt5a', 'Yjefn3', 'Zfpm2'};

% function GeneList = FavoriteGenes();
% GeneList = {'Adcy1', 'Adarb2', 'Adra1b', 'Alcam', 'Arl4c', 'Arx', 'Asic4', 'Atp6v1g1', 'Bmper', ...
%     'Cacna1a', 'Cacna1b', 'Cacna1e', 'Cacna2d1', 'Cacna2d2', 'Cacna2d3', 'Cacnb2', 'Cacng2', 'Cadps2', 'Calb1', 'Calb2', 'Caln1', 'Car2', 'Cck', 'Cckbr', 'Cd24a', 'Cdh13', 'Cend1', ...
%     'Chrm1', 'Chrm2', 'Chrm3', 'Chrna4', 'Clstn2', ...
%     'Cnr1', 'Cntnap2', 'Col19a1', 'Col25a1', 'Cplx1', 'Cplx2', 'Cplx3', 'Cpne5', 'Crabp1', 'Crh', 'Crhbp', 'Cryab', 'Csrp2', 'Cttnbp2', ...
%     'Cxcl14', 'Dlx1os', 'Dlx6os1', 'Dpy19l1', 'Elfn1', 'Enpp2', 'Epha4', 'Fxyd6', 'Gabra1', 'Gabra3', 'Gabrd', ...
%     'Gad1', 'Gad2', 'Gda', 'Gnb4', 'Gng2', 'Gria3', 'Grik1', 'Grik5', 'Grin3a', 'Grip1', 'Grm1', 'Gstp1', 'Gucy1a3', ...
%     'Hapln1', 'Hpcal1', 'Htr2c', 'Htr3a', 'Id2', 'Igf1', 'Igfbp5', 'Itm2b', 'Kcna1', 'Kcna2', 'Kcnc1', 'Kcnc2', 'Kcnh1', 'Kcnip1', 'Kcnip2', 'Kcnmb2', 'Kcnmb4',  'Kctd12', 'Kirrel3', 'Kit', 'Krt73', ...
%     'Lamp5', 'Ldb2', 'Lgals1', 'Lhx6', 'Lin7a', ...
%     'Lypd6','Lypd6b', 'Mafb', 'Mef2c', 'Mpped1', 'Ncald', 'Ndnf', 'Necab1', 'Nell1', ...
%     'Nfib', 'Npy', 'Npas1', 'Npas3', 'Npy1r', 'Npy2r', 'Nr2f2',  'Nrsn1', 'Nrn1', 'Nsg1', 'Ntng1', 'Nxph1', ...
%     'Olfm3', 'Parm1', 'Pbx1', 'Pcdh8', 'Pcp4', 'Penk', 'Phf20l1', 'Plcxd3', 'Pls3', 'Pnoc', 'Ppargc1a', 'Prkca', 'Prkacb', 'Psd3', 'Ptprd', 'Ptpre', 'Ptrhd1', 'Pvalb', ...
%     'Rab3a', 'Rab3b', 'Rbfox1', 'Rbp4', 'Reln', 'Resp18', 'Rgs8', 'Rgs10', 'Rgs12', 'Rgs16', 'Rgs17', 'Ryr3', ...
%     'Satb1', 'Scn1a', 'Scn8a', 'Sec61b', 'Sel1l3', 'Sema3c', 'Sema5a', 'Serpine2','Sez6',  ...
%     'Slc16a2', 'Slc17a8', 'Slc32a1', 'Slc38a1', 'Slc6a1', 'Slc7a14', 'Snca', ...
%     'Sox5', 'Sox6', 'Sparcl1', 'Spag9', 'Spock1', 'Srgap3', 'Sst', 'Sstr1', 'Stk32a', 'Strip2', 'Syt6', 'Tac1', 'Tac2', 'Tacr1', 'Tenm1', 'Tgfb2', 'Thrsp', ...
%     'Tnfaip8l3', 'Trp53i11', 'Trhde', 'Trpc6', 'Vamp1', 'Vip', 'Vwa5a', 'Wnt5a', 'Zbtb18', 'Zbtb20', 'Zfpm2'};

return
% This has more
% GeneList = {'Adcy1', 'Adarb2', 'Adra1b', 'Alcam', 'Aldoa', 'App', 'Arl4c', 'Arx', 'Asic4', 'Atp6v1g1', 'Atp1b1', 'Atp2b1', 'B3galt2', 'Bmper', ...
%     'Cacna1a', 'Cacna1b', 'Cacna1e', 'Cacna2d1', 'Cacna2d2', 'Cacna2d3', 'Cacnb2', 'Cacng2', 'Cadps2', 'Calb1', 'Calb2', 'Caln1', 'Camk2a', 'Camk2b', 'Car2', 'Cck', 'Cckbr', 'Cd200', 'Cd24a', 'Cdh13', 'Cend1', ...
%     'Chrm1', 'Chrm2', 'Chrm3', 'Chrna4', 'Clstn2', ...
%     'Cnr1', 'Cntnap2', 'Col19a1', 'Col25a1', 'Cplx1', 'Cplx2', 'Cplx3', 'Cpne5', 'Crabp1', 'Crh', 'Crhbp', 'Cryab', 'Csrp2', 'Cttnbp2', ...
%     'Cxcl14', 'Dlx1os', 'Dlx6os1', 'Dpy19l1', 'Elfn1', 'Enpp2', 'Epha4', 'Fxyd6', 'Gabra1', 'Gabra3', 'Gabrd', ...
%     'Gad1', 'Gad2', 'Gda', 'Gnb4', 'Gng2', 'Gria3', 'Grik1', 'Grik5', 'Grin3a', 'Grip1', 'Grm1', 'Grm7', 'Gstp1', 'Gucy1a3', ...
%     'Hapln1', 'Hpcal1', 'Htr2c', 'Htr3a', 'Id2', 'Igf1', 'Igfbp5', 'Itm2b', 'Kcna1', 'Kcna2', 'Kcnc1', 'Kcnc2', 'Kcnh1', 'Kcnip1', 'Kcnip2', 'Kcnmb2', 'Kcnmb4',  'Kctd12', 'Kirrel3', 'Kit', 'Krt73', ...
%     'Lamp5', 'Ldb2', 'Lgals1', 'Lhx6', 'Lin7a', 'Luc7l3', ...
%     'Lypd6','Lypd6b', 'Mafb', 'Mcu', 'Mef2c', 'Mpped1', 'Ncald', 'Ndnf', 'Necab1', 'Nell1', ...
%     'Nfib', 'Npy', 'Npas1', 'Npas3', 'Npy1r', 'Npy2r', 'Nr2f2',  'Nrip1', 'Nrsn1', 'Nrn1', 'Nsg1', 'Ntng1', 'Nxph1', ...
%     'Olfm3', 'Parm1', 'Pbx1', 'Pcdh8', 'Pclo', 'Pcp4', 'Penk', 'Phf20l1', 'Plcxd3', 'Pls3', 'Pnoc', 'Ppargc1a', 'Prkca', 'Prkacb', 'Psd3', 'Ptprd', 'Ptpre', 'Ptrhd1', 'Pvalb', ...
%     'Rab3a', 'Rab3b', 'Rbfox1', 'Rbp4', 'Reln', 'Resp18', 'Rgs8', 'Rgs10', 'Rgs12', 'Rgs16', 'Rgs17', 'Ryr3', ...
%     'Satb1', 'Scn1a', 'Scn8a', 'Sec61b', 'Sel1l3', 'Sema3c', 'Sema5a', 'Serpine2','Sez6',  ...
%     'Slc12a2', 'Slc16a2', 'Slc17a8', 'Slc32a1', 'Slc38a1', 'Slc6a1', 'Slc7a14', 'Snca', ...
%     'Sox5', 'Sox6', 'Sparcl1', 'Spag9', 'Spock1', 'Srgap3', 'Sst', 'Sstr1', 'Stk32a', 'Strip2', 'Syt6', 'Tac1', 'Tac2', 'Tacr1', 'Tenm1', 'Tgfb2', 'Thrsp', 'Tm2d1', ...
%     'Tnfaip8l3', 'Trp53i11', 'Trhde', 'Trpc6', 'Usp28', 'Vamp1', 'Vip', 'Vwa5a', 'Wnt5a', 'Zbtb18', 'Zbtb20', 'Zfpm2', 'Zwint'};

%return

    'Cnr1'
    'Vip'
    'Cadps2'
    'Kctd12'
    'Cck'
    'Rgs12'
    'Sncg'
    'Trp53i11'
    'Fxyd6'
    'Cxcl14'
    'Htr3a'
    'Cplx2'
    'Sln'
    'Necab1'
    'Serpine2'
    'Col25a1'
    'Krt73'
    'Npas1'
    'Adarb2'
    'Cpne5'
    'Gng2'
    'Igf1'
    'Snca'
    'Cst3'
    'Necab2'
    'Fabp5'
    'Col19a1'
    'Tac2'
    'Eps8'
    'Igfbp5'
    'Syt6'
    'Car2'
    'Elmo1'
    'Ubash3b'
    'Sst'
    'Crhbp'
    'Hapln1'
    'Lamp5'
    'Sparcl1'
    'Lgals1'
    'Nxph1'
    'Nrsn1'
    'Ncald'
    'Satb1'
    'Grin3a'
    'Adcy1'
    'Rgs17'
    'Lhx6'
    'Npy'
    'Prkacb'
        'Vip'
    'Tac2'
    'Krt73'
    'Kctd12'
    'Atp2b1'
    'Alcam'
    'Rgs10'
    'Prkcb'
    'Nrgn'
    'Epha4'
    'Sparcl1'
    'Slc17a8'
    'Crym'
    'Enpp2'
    'Npy'
    'Sema3c'
    'Cd24a'
    'Nrip3'
    'Tmsb4x'
    'Rgs8'
    'Nov'
    'Reln'
    'Igfbp5'
    'Ndnf'
    'Spock1'
    'Rgs12'
    'Cpne5'
    'Cxcl14'
    'Serpine2'
    'Vsnl1'

    'Crym'
    'Rgs4'
    'Hpca'
    'Ppp3ca'
    'Sv2b'
    'Prkcb'
    'Gda'
    'Epha4'
    'Pcdh17'
    'Kctd12'
    'Npy'
    'Krt73'
    'Rgs10'
    'Cck'
    'Slc17a8'
    'Car2'
    'Sncg'
    'Tac2'
    'Dlx1os'
    'Cnr1'
    'Slc2a13'
    'Trp53i11'
    'Id2'
    'Elmo1'
    'Mgll'
    'Glul'
    'Sema3c'
    'Fabp5'
    'Enpp2'
    'Abat'
    'Htr3a'
    'Ubash3b'
    'Slc22a23'
    'Cd24a'
    'Necab2'
    'Rasgrf2'
    'Lars2'
    'Pclo'
    'Cadps2'
    'Arl4c'
    'Alcam'
    'Col25a1'
    'Dlx6os1'
    'Maf'
    'Shisa9'
    'Slc6a1'
    'Sln'
    'Cntnap2'
    'Osbpl8'
    'Stk32c'
    'Elavl2'
    'Kcnip1'
    