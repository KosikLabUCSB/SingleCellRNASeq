BrainTREM2 <- seurat_combined

library(Seurat)


#make list of celltype markers
Ast<- c("AQP4","GJA1","GJB6","SLC4A4","SLC1A2","F3","BMPR1B","FGFR3",
        "SLC39A12","CLDN10","DIO2","ALDOC","ALDH1L1","SLC1A3","CLU",
        "ATP13A4","SLCO1C1","SLC14A1","CHRDL1","GPR37L1","ACSBG1",
        "ATP1A2","SLC25A18","EDNRB","PPAP2B","GFAP","SOX9","SDC4",
        "PPP1R3C","NCAN","MLC1","GLI3","SLC7A11","ACSL6","RFX4","ID4",
        "AGT","SFXN5","GABRG1","PAX6","RORB","GRM3","PTPRZ1","PSD2",
        "SLC6A11","ATP1B2","NTSR2","S1PR1","SLC15A2","ELOVL2","TRIL",
        "SCARA3","MGST1","KIAA1161","FAM107A","BCAN","SPARCL1","NWD1",
        "NTRK2","SLC7A10","SCG3","ACOT11","KCNN3","MFGE8","RANBP3L",
        "GPC5","EZR","ADHFE1","GABRB1","TMEM47","PAMR1","CPE","FABP7",
        "LIX1","SLC13A5","IL33","SLC7A2","EGFR","PREX2","NDRG2","DTNA",
        "ABCD2","HEPACAM","RGS20","ARHGEF26","GPAM","CHI3L1","ADCYAP1R1",
        "GDPD2","SLC1A4","POU3F2","ETNPPL","MEGF10","MT3","TTYH1",
        "PRODH","PLCD4","DDAH1","LGR4","HTRA1")

Endo<- c("APOLD1","FLT1","RGS5","PTPRB","TM4SF1","ABCB1","ITM2A","SDPR",
         "SLCO1A2","FN1","EMCN","ESAM","NOSTRIN","CD34","SLC38A5","CYYR1",
         "PODXL","CDH5","VWF","MECOM","CD93","ABCG2","TEK","PALMD","ERG",
         "CLDN5","PECAM1","KDR","ITGA1","ICAM2","ATP10A","ANXA3","CA4",
         "MYCT1","GIMAP6","ANXA1","PTRF","KIAA1462","EBF1","HMCN1","ENG",
         "IGFBP7","ARHGAP29","ANXA2","OCLN","HIGD1B","SLC2A1","GNG11",
         "SLC19A3","EPAS1","TBX3","SRGN","SOX7","SLC16A4","CAV1","CLIC5",
         "VIM","HEG1","CCDC141","C10ORF10","EDN1","ROBO4","TMEM204",
         "PROM1","IFITM1","LEF1","COBLL1","WWTR1","HBB","ETS1","SLC39A8",
         "COL4A1","OSMR","ADCY4","TIE1","EDN3","THBD","BSG","AHNAK",
         "MYO1B","IL1R1","CXCL12","CLEC14A","GATA2","SGPP2","SHE","PLTP",
         "SPARC","ACVRL1","MMRN2","NID1","TNFSF10","FOXC1","UACA","CGNL1",
         "MFSD2A","NET1","ABCC9","FLI1","C1ORF54")

Mic<- c("CCL4","CCL3","CSF1R","CX3CR1","P2RY12","C1QB","RGS1","GPR183",
        "GPR34","CTSS","LAPTM5","CD53","IL1A","C3AR1","PLEK","FCGR2A",
        "CD83","ITGAM","P2RY13","CD86","TREM2","TYROBP","FCER1G","NCKAP1L",
        "SELPLG","SLC2A5","CD14","C1QC","C1QA","MPEG1","HAVCR2","PTAFR",
        "LY86","AIF1","ALOX5AP","LPCAT2","SLA","PTPRC","FCGR1A","CCL2",
        "BLNK","IL10RA","BCL2A1","C5AR1","RHOH","CD84","CSF3R","TLR7",
        "TLR2","HPGDS","LCP1","CD300A","FYB","MRC1","FAM105A","IRF8",
        "LCP2","RGS10","CD74","PTPN6","TBXAS1","LYZ","DOCK2","TMEM119",
        "NLRP3","ARHGDIB","CCRL2","IKZF1","ARHGAP25","DOCK8","HEXB",
        "THEMIS2","SAMSN1","HK2","PLD4","APBB1IP","ITGB2","RUNX1",
        "SLCO2B1","TLR1","FGD2","HCLS1","GPR84","OLFML3","MAFB","PIK3CG",
        "SIGLEC7","IL1B","PIK3R5","IL6R","CXCL16","CLEC4A","PTGS1","SUSD3",
        "LYN","VAV1","SLC11A1","RBM47","SYK","C10ORF128")

Neu<- c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3",
        "GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST",
        "VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5",
        "NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK",
        "ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A",
        "RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87",
        "ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2",
        "DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2",
        "CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4",
        "ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2",
        "CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1",
        "GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2",
        "RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")


Oli<- c("PLP1","MOBP","CLDN11","MBP","UGT8","ERMN","MOG","MAG","OPALIN","CNP",
        "MAL","GPR37","TF","MYRF","GJB1","ASPA","ENPP2","BCAS1","LPAR1","FA2H",
        "ENPP6","APOD","CNTN2","CRYAB","KLK6","ERBB3","ANLN","SEPT4","PLEKHB1",
        "TMEFF2","ST18","PTGDS","PEX5L","SLAIN1","QDPR","PLLP","TMEM125","HHIP",
        "LGI3","TUBB4A","PLEKHH1","S1PR5","MAP6D1","GSN","EVI2A","EDIL3",
        "CMTM5","GJC3","CA14","NFASC","TPPP","TMEM88B","TRIM59","CDH19","APLP1",
        "NIPAL4","ADAMTS4","STMN4","S100B","CA2","PRR18","OLIG1","FOLH1",
        "NINJ2","NDRG1","SLC24A2","SGK2","GALNT6","KCNA1","SH3TC2","TTLL7",
        "SH3GL3","DOCK5","SCD","FEZ1","SLC44A1","RHOU","PPP1R16B","TSPAN2",
        "C10ORF90","TNFAIP6","NKAIN2","MOB3B","PRKCQ","PPP1R14A","PLA2G16",
        "DBNDD2","CDK18","PCDH9","ANO4","AGPAT4","OMG","FGFR2","TMEM63A",
        "GLTP","CCP110","PLEKHG3","RAB33A","PSAT1","ZNF536")

OPCs<-c("PDGFRA","TNR","PCDH15","SHC4","VCAN","LHFPL3","NEU4","GPR17",
        "PTPRZ1","OLIG1","MMP16","DSCAM","C8ORF46","SEMA5A","MATN4",
        "UGT8","GRIA3","CNTN1","BCAS1","SULF2","LUZP2","GJC3","NXPH1",
        "APOD","MEGF11","LRRTM3","BRINP3","GALNT13","GRIA4","MYT1","SUSD5",
        "LRRN1","SOX10","PRKCQ","SOX6","ITGB8","TMEM255A","GFRA1","RLBP1",
        "PNLIP","XYLT1","GPSM2","TMEM255B","SEZ6L","STK32A","C14ORF37",
        "LPPR5","SEMA3D","CSPG4","CSMD3","TMEM132B","SCRG1","KCNH8",
        "CACNG4","UGDH","DPP6","BCAT1","PLLP","ERBB3","RNF43","S100B",
        "SORCS1","OLIG2","CHRNA4","KCNJ16","PPAPDC1A","CSMD1","OPCML",
        "PRKG2","COBL","FIGN","ACAN","TGFA","NLGN1","SLC6A13","EMID1",
        "CHST6","TMEM100","GAL3ST1","EDIL3","KCNJ10","SLITRK3","SNTG1",
        "CSPG5","ERBB4","SLC35F1","B3GAT2","C1QL1","SERINC5","CKAP2",
        "LRRTM4","DPYD","SLITRK1","NCALD","CALCRL","SPP1","ZNF488",
        "ADAM12","SULF1","HAS2")

Peri<-c("KDR","APLNR","MFAP4","MCAM","DDR2","COL5A2","MYO1B","TNS1","LMOD1",
        "COL3A1","IGFBP5","DES","FMOD","PRRX1","DPT","RGS5","HSD11B1","VIM",
        "FAP","FBN1","ANGPTL2","SLC38A11","FRZB","SERPING1","NR1H3","PAMR1",
        "FIBIN","ADAM33","REM1","COX4I2","GNB4","PDE5A","SYNPO2","FABP4",
        "PCDH18","POSTN","P2RY14","ECM1","OLFML3","COL15A1","SVEP1","TEK",
        "COL16A1","PODN","HEYL","ALPL","HSPB7","MXRA8","LDB2","STEAP4","HTRA3",
        "MSX1","SOD3","PDGFRA","SPARCL1","ELN","ABCC9","COL1A2","RARRES2",
        "AQP1","FBLN2","CXCL12","C1S","KCNJ8","COG7","COX7A1","HSPB6","MFGE8",
        "ANPEP","ADM","HTRA1","IFITM1","MRGPRF","COL4A2","COL4A1","EDNRA",
        "INPP4B","CDH11","ANGPT2","MMP2","CDH5","TAGLN","PTH1R","ICAM1","ROBO4",
        "THY1","C1QTNF5","HSPB2","CRYAB","CSPG4","ISLR","RASL12","TBX18","NT5E",
        "ZIC1","PLSCR4","BGN","LAMA2","LAMA4","COL6A2","COL6A1","TIMP3","DCN",
        "LUM","GLI1","NDUFA4L2","SGCD","SLIT3","ADAMTS2","PECAM1","EMID1",
        "SPARC","VTN","COL1A1","RAMP2","AOC3","HIGD1B","CYGB","C1QTNF1",
        "COLEC11","CLEC14A","DLK1","FOXC1","ASPN","ECM2","PDZD2","ANGPT1",
        "COL14A1","GPIHBP1","ACVRL1","ATP13A5","FSTL1","ADAMTS5","FNDC1",
        "THBS2","DAAM2","NOTCH3","VEGFA","PDGFRB","CD248","EFEMP2","ACTA2",
        "COL6A3","ABCA8","MXRA5","C1R","HSPB2-C11ORF52","MEDAG","TNS2","ADGRA2",
        "ADGRL4","ADGRF5","CCN2","PLPP3","LHFPL6")

Epen<- c("PLXNB2","PLTP","C20ORF85","EFHC1","APOBEC4","PCP4L1","PTPRT","ARMC3",
         "ENKUR","MORN5","TTLL9","MEIG1","ANGPTL2","WDR38","HDC","SPEF1",
         "TMEM212","ZBBX","SPAG17","ADH7","STOML3","TM4SF1","CELSR2","FAM166B",
         "TEKT2","CROCC","CHEK2","CCDC60","CRYGN","RARRES2","AQP1","USP18",
         "LRRC23","VWA3A","CXCL17","MYO16","TTC29","NWD1","HYDIN","TMEM231",
         "DYNLRB2","TRIM71","CCDC153","CALML4","ZMYND10","AKAP14","LRRIQ1",
         "MYB","RSPH4A","PPIL6","S100B","TMEM107","EFNB3","TEKT1","GAS2L2",
         "KRT15","GFAP","FOXJ1","CCDC40","AK7","HSPA2","SNTN","CAPSL","SPEF2",
         "PVALB","ODF3B","SPAG6","IQCG","EFCAB1","PACRG","TMEM232","TEKT4",
         "RSPH1","SPATS1","TCTE1","RFX2","SIX3","NME5","AQP4","SCGB1A1",
         "LRRC10B","RIIAD1","CYP2A13","DNAH7","DNAH9","C1ORF87","DNAI1",
         "C2ORF73","RGS22","DNAH12","FAM183A","WFDC6","DNAI2","C1ORF194",
         "DNAH6","BPIFA1","C5ORF49","DNAAF1","DNAH2","FAM81B","MAP3K19",
         "VWA3B","AK8","C9ORF24","C6ORF118","LRRC71","PPP1R42","NME9","BPIFB1",
         "MS4A8","DRC1","PIFO","EPPIN-WFDC6","CATIP","TEX26","FAM216B",
         "C11ORF97","ANKRD66","CFAP44","CFAP70","ERICH3","CFAP43","CFAP126",
         "CFAP61","CFAP54","CFAP221","CFAP52","CFAP99","DRC7","CFAP47","SAXO2",
         "CFAP100","CFAP77","CFAP65","DRC3","CFAP206")




DefaultAssay (BrainTREM2) <- "RNA"

#BrainSUBSET <- subset(BrainFP.combined, features=c(Ast, Endo, Mic, Neu, Oli, OPCs, Peri))

celltype<-'Peri'

BrainSUBSET <- subset(BrainTREM2, features=get(celltype))

av.expBrainSUBSET <- AverageExpression(BrainSUBSET)$SCT

cal_z_score <- function(x){(x - mean(x)) / sd(x)}

av.expBrainSUBSET <- t(apply(av.expBrainSUBSET, 1, cal_z_score))

library(pheatmap)
plot<-pheatmap (av.expBrainSUBSET, cluster_rows = F, fontsize_row = 7, main= paste(celltype))



#Rename the cluster----
###Rename all cluster names (replace the names for whatever you have identified)
new.cluster.ids<- c('Oli','Oli','Oli','Oli','Oli','Neu.2','Mic','Ast','Oli','Oli','Ast','Neu.2','Neu.2','Neu.2','Oli','Neu.1','Neu.1','Mic','Neu.2','Neu.2','Oli','Neu.2','OPC','OPC','End','Ast','Ast','Neu.2','Neu.2','End/Per','Mic','Neu.2','Ast','Oli','Neu.2','Neu.2','Neu.2','Neu.2','End/Per','Unk','Unk')
  
  
names(new.cluster.ids) <- levels(BrainTREM2)
BrainTREM2 <- RenameIdents(BrainTREM2, new.cluster.ids)

DimPlot(BrainTREM2, label = F) 

BrainTREM2[["Heatmap.celltype"]] <- Idents(object = BrainTREM2)

saveRDS(BrainTREM2, "~/TREM2/BrainTREM2.rds")


#remove MT genes from counts matrix

##############################################################################################################################################################################
#### Remove all MT genes
####Get the total genes expressed list from seurat object
total.genes <- list(rownames(BrainTREM2))
total.genes <- as.data.frame(do.call(cbind, total.genes))
mito.genes <- grep(pattern = "^MT-", x = rownames(x = BrainTREM2), value = TRUE)

####get the MT-genes rows in the matrix
MTGenes <- which(toupper(rownames(BrainTREM2)) %in% mito.genes)
MTGenes #### check which rows your MT genes are

total.genes_subset <- total.genes[-c(MTGenes),] ###add the rows according to the print of your MTGenes function
total.genes_subset <- as.character (total.genes_subset) 

### Subset Seurat object containg all gene reads but the MT genes 
BrainTREM2 <- subset (BrainTREM2, features = total.genes_subset)
#######################################################################################
saveRDS(BrainTREM2, "~/TREM2/BrainTREM2-noMT.rds")
