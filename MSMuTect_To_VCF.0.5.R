cat(paste0(substr(Sys.time(), 1, 19), "  @@@ MSMuTect0.5 to VCF - start:\n"))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# built-in filters:
# variants in allosomes (X&Y) are removed.
# Duplicated variants are removed
# variants in which the tumor allele is identical to the reference genome are removed.

#selected_arguments:
cat(paste0("\tRaw input:",args,"\n"))
cat("\n")
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ get arguments:\n"))
#TODO set default values
if (length(args) != 7) {
  stop(
paste0("\n7 argument must be supplied: (you entered",length(args),")\n
1-input file (char)\n
2-InputPath (char)\n
3-Output Path (char)\n
4-Single Sample Mode (Logical T/F)\n
5-Filtering mode (int 1..7):\n
\t1. Filter all non REF normals,\n
\t2. Filter normals homozygous to ALT,\n
\t3. Filter normals heterozygous to ALT,\n
\t4. No filters,\n
\t5. Exclude LOH: Remove variants which the normal is heterozygous and the tumor homozygous to one of the normal alleles.\n
\t6. Remove all normals with heterozygous samples.\n
\t7. Remove all normals with heterozygous or homozygous ALT alleles & all tumors that are not diploid heterozygous to REF.\n
\t8. keep only tumor alleles that contain novel variants(not seen in REF or Normal./n)
6-output-filtered MSMuTect file together with your vcf (Logical T/F)\n
7-Split multiAllelic loci in .vcf to single lines (Logical T/F)\n\n"),
    call. = FALSE
  )
}
# example input args<-c("LP5007.A-26-A1_T._vs_.A-26-A1_N.full.mut.MutOnly","C:/Users/Bad-G/OneDrive - Technion/Third_Degree/MSMuTect_To_VCF_v0.5","C:/Users/Bad-G/Desktop","F","8",F,F)
INPUT <- args[1]
inputPath <- args[2]
outputPath <- args[3]
OneSample <- as.logical(args[4])
mode <- as.integer(args[5])
output_MSMuTect <- as.logical(args[6])
singleLine <- as.logical(args[7])
cat("Processed input:\n")
cat(paste0("\t",1,") ","Input file:.......................\t", INPUT,"\n"))
cat(paste0("\t",2,") ","Input Path:.......................\t", inputPath,"\n"))
cat(paste0("\t",3,") ","OutputPath:.......................\t", outputPath,"\n"))
cat(paste0("\t",4,") ","singleSample Mode:................\t", OneSample,"\n"))
cat(paste0("\t",5,") ","Normal filter mode:...............\t", mode,"\n"))
cat(paste0("\t",6,") ","Output MSMuTect file with filters:\t", output_MSMuTect,"\n"))
cat(paste0("\t",7,") ","Separate multiallelic variants:...\t",singleLine,"\n\n"))
#additional arguments
remove_Sex_Chr <- T

#different modes to filter variants in normal sample
MODES <- c(
  "1) Filter all non REF normals ",
  "2) Filter normals homozygous to ALT ",
  "3) Filter normals heterozygous to ALT ",
  "4) No filters ",
  "5) Exclude LOH: Remove variants which the normal is heterozygous and the tumor homozygous to one of the normal alleles. ",
  "6) Remove all normals with heterozygous samples. ",
  "7) Remove all normals with heterozygous or homozygous ALT alleles & all tumors that are not diploid heterozygous to REF. ",
  "8) Keep only tumor alleles that contain novel variants (not seen in REF or Normal). "
)
if (output_MSMuTect & !singleLine){
  warning("Cant export MSMuTect file when split multiallelic option is on (argument 7==T).")
}
### --- import MSMuTect input  --- ###
setwd(inputPath)
options(stringsAsFactors = FALSE)
MSMuTect_Input <- read.delim(INPUT)
MSMuTect_output<-MSMuTect_Input[!is.na(MSMuTect_Input$TUMOR_ALLELE_1),] #Fix Error in MSMuTect0.5 some calls have NA in their tumor!
rm(MSMuTect_Input)
invisible(gc())
### --- generate .vcf file --- ###
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ generate VCF header\n"))

VCF_output<-as.data.frame(matrix(NA,ncol = 9,nrow = nrow(MSMuTect_output) ))
colnames(VCF_output)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
if(OneSample){
  VCF_header<-paste0(
    '##fileformat=VCFv4.2
##fileDate=',Sys.Date(),'
##INFO=<ID=LOCUS,Number=.,Type=String,Description="Microsatelite locus taken from MSMuTect output">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##source=Yosef Maruvka`s lab:www.MaruvkaLab.com
##Tool=github.com/Hagay-Ladany/MSMuTect-0.5-to-VCF Input:',INPUT,',',inputPath,',',outputPath,',',OneSample,',',mode,',',output_MSMuTect,',',singleLine,'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR')
  
}else{
  VCF_header<-paste0(
    '##fileformat=VCFv4.2
##fileDate=',Sys.Date(),'
##INFO=<ID=LOCUS,Number=.,Type=String,Description="Microsatelite locus taken from MSMuTect output">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##Tool=github.com/Hagay-Ladany/MSMuTect-0.5-to-VCF Input:',INPUT,',',inputPath,',',outputPath,',',OneSample,',',mode,',',output_MSMuTect,',',singleLine,'
##source=Yosef Maruvka`s lab:www.MaruvkaLab.com
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR')
}

cat(paste0(substr(Sys.time(),1,19),"  @@@ generate VCF columns\n"))

VCF_output[, "#CHROM"] <- MSMuTect_output$CHROMOSOME
VCF_output[, "POS"]    <- MSMuTect_output$START
VCF_output[, "ID"]     <- "."
VCF_output[, "QUAL"]   <- "."
VCF_output[, "FILTER"] <- "PASS"
VCF_output[,"INFO"]    <- paste0("LOCUS=",with(MSMuTect_output,paste(CHROMOSOME,START,END,PATTERN,floor(REFERENCE_REPEATS),sep =":")))
VCF_output[, "FORMAT"] <- "GT"

REFERENCE_REPEATS<-floor(MSMuTect_output$REFERENCE_REPEATS)

# generate REF and ALT columns
cat(paste0(substr(Sys.time(),1,19),"  @@@ generate VCF ALT REF columns\n"))

if (OneSample){ #check if MSMuTect  run was on single sample mode
  GT<-as.matrix(cbind(REFERENCE_REPEATS,MSMuTect_output[,c("ALLELE_1","ALLELES_2","ALLELES_3","ALLELES_4")]))
  for (i in 2:5){GT[GT[,i]==GT[,1],i]<-NA} #remove REF identical
  RemoveMe<-!sapply(1:nrow(GT),function(X) all(is.na(GT[X,2:5])) )
  GT<-GT[RemoveMe,] # remove rows where all alleles are 'NA'
  
  #fix of tables:
  VCF_output<-VCF_output[RemoveMe,]
  MSMuTect_output<-MSMuTect_output[RemoveMe,]
  REFERENCE_REPEATS<-REFERENCE_REPEATS[RemoveMe]
  
  GT<-as.data.frame(GT)
  GT<-t(apply(GT,1,function(X) X - (min(X,na.rm = T)-1))) # shrink repeats
  if(nrow(GT)>1){
    GT_unique<-apply(GT[,2:5],1,function(x) unique(na.omit(x)))
  }else{
    GT_unique<-unique(na.omit(GT[,2:5]))
  }
  GT<-cbind(GT[,"REFERENCE_REPEATS"],substr(MSMuTect_output$REFERENCE_SEQUENCE,1,nchar(MSMuTect_output$PATTERN)))#get ref patterns
  colnames(GT)<-c("REF","Base")
  ALT_ALLELS<-sapply(1:length(GT_unique),function(ROW) strrep(GT[ROW,"Base"], GT_unique[[ROW]]))
  VCF_output[,"ALT"]<-sapply(ALT_ALLELS,function(x) paste0(na.omit(x),collapse = ","))
  # MSMuTect motif detection may by different from what's really found in the reference genome used in the following steps - in case MSMuTect Reference genome equal to the requested genome we can use the ref sequence in out mutect file
  VCF_output[,"REF"]<-sapply(1:nrow(GT),function(ROW) substr(MSMuTect_output$REFERENCE_SEQUENCE[ROW],1,nchar(GT[ROW,"Base"])*as.numeric(GT[ROW,"REF"]))  )
  #table((nchar(VCF_output$REF)-nchar(VCF_output$ALT))>0)  check del/dup ratio
}else{ #regular run (Normal / Tumor)
  GT<-as.matrix(cbind(REFERENCE_REPEATS,MSMuTect_output[,c("NORMAL_ALLELE_1","NORMAL_ALLELES_2","NORMAL_ALLELES_3","NORMAL_ALLELES_4","TUMOR_ALLELE_1","TUMOR_ALLELES_2","TUMOR_ALLELES_3","TUMOR_ALLELES_4")]))
  for (i in 2:9){GT[GT[,i]==GT[,1],i]<-NA} #remove REF identical
  
  GT<-as.data.frame(GT)
  GT<-t(apply(GT,1,function(X) X - (min(X,na.rm = T)-1))) # shrink repeats
  if(nrow(GT)>1){
    GT_unique<-apply(GT[,2:9],1,function(x) unique(na.omit(x)))
  }else{
    GT_unique<-unique(na.omit(GT[,2:9]))
  }
  GT<-cbind(GT[,"REFERENCE_REPEATS"],substr(MSMuTect_output$REFERENCE_SEQUENCE,1,nchar(MSMuTect_output$PATTERN))) ; colnames(GT)<-c("REF","Base")#cut motif from normal allele (so different ALT could be created)
  ALT_ALLELS<-sapply(1:length(GT_unique),function(ROW) strrep(GT[ROW,"Base"], GT_unique[[ROW]]))
  VCF_output[,"ALT"]<-sapply(ALT_ALLELS,function(x) paste0(na.omit(x),collapse = ","))
  # MSMuTect motif detection may by different from what's really found in the reference genome used in the following steps - in case MSMuTect Reference genome equal to the requested genome we can use the ref sequence in out mutect file
  VCF_output[,"REF"]<-sapply(1:nrow(GT),function(ROW) substr(MSMuTect_output$REFERENCE_SEQUENCE[ROW],1,nchar(GT[ROW,"Base"])*as.numeric(GT[ROW,"REF"]))  )
}
#give score for each sample in N/N format
#write ALT to VCF
if (OneSample){
  GT<-as.matrix(cbind(REFERENCE_REPEATS,MSMuTect_output[,c("ALLELE_1","ALLELES_2","ALLELES_3","ALLELES_4")]));
  GT_unique<-apply(GT[,1:5],1,function(x) unique(na.omit(x)));
  GT<-GT[,-1];
  GT<-as.data.frame(GT);
  GT_score<-t(sapply(1:nrow(GT),function(x) match(GT[x,] ,GT_unique[[x]])-1));
  colnames(GT_score)<-colnames(GT);
  #assuming diploid chrs - (will be wrong for sex chromosomes!) I duplicate the first allele to the missing second ones to achive at least diploid chrs (even though biologically this is not always true)
  GT_score[is.na(GT_score[,2]),2]<-GT_score[is.na(GT_score[,2]),1]; #fill-up empty second allele with the first allele content.
  GT_Normal<-paste(GT_score[,1],GT_score[,2],GT_score[,3],GT_score[,4],sep = "/");
  GT_Normal<-gsub("/NA","",GT_Normal); #remove NA's from tumor genotypes
  VCF_output[,"TUMOR"]<-GT_Normal
  
}else{
  GT<-as.matrix(cbind(REFERENCE_REPEATS,MSMuTect_output[,c("NORMAL_ALLELE_1","NORMAL_ALLELES_2","NORMAL_ALLELES_3","NORMAL_ALLELES_4","TUMOR_ALLELE_1","TUMOR_ALLELES_2","TUMOR_ALLELES_3","TUMOR_ALLELES_4")]));
  GT_unique<-apply(GT[,1:9],1,function(x) unique(na.omit(x)));
  GT<-GT[,-1];
  GT<-as.data.frame(GT);
  GT_score<-t(sapply(1:nrow(GT),function(x) match(GT[x,] ,GT_unique[[x]])-1));
  colnames(GT_score)<-colnames(GT);
  #assuming diploid chrs - (will be wrong for sex chromosomes!) I duplicate the first allele to the missing second ones to achive at least diploid chrs (even though biologically this is not always true)
  GT_score[is.na(GT_score[,2]),2]<-GT_score[is.na(GT_score[,2]),1]; #fill-up empty second allele with the first allele content.
  GT_Normal<-paste(GT_score[,1],GT_score[,2],GT_score[,3],GT_score[,4],sep = "/");
  GT_Normal<-gsub("/NA","",GT_Normal); #remove NA's from tumor genotypes
  GT_score[is.na(GT_score[,6]),6]<-GT_score[is.na(GT_score[,6]),5]; #fill-up empty second allele with the first allele content.
  GT_Tumor <-paste(GT_score[,5],GT_score[,6],GT_score[,7],GT_score[,8],sep = "/");
  GT_Tumor <-gsub("/NA","",GT_Tumor) #remove NA's from tumor genotypes
  GT_Tumor[nchar(GT_Tumor)==1]<-paste0(GT_Tumor[nchar(GT_Tumor)==1],"/",GT_Tumor[nchar(GT_Tumor)==1])  #assuming minimum of diploid chrs in tumor as well.
  VCF_output[,"NORMAL"]<-GT_Normal
  VCF_output[,"TUMOR"] <-GT_Tumor
}
#remove sex chrs?
if (remove_Sex_Chr) {
  if (output_MSMuTect & !singleLine){
    MSMuTect_output <-
      MSMuTect_output[VCF_output[, "#CHROM"] != "X" &
                        VCF_output[, "#CHROM"] != "Y", ]
  }
  VCF_output <-
    VCF_output[VCF_output[, "#CHROM"] != "X" &
                 VCF_output[, "#CHROM"] != "Y", ]
}

#order
if (output_MSMuTect & !singleLine){
  MSMuTect_output <-MSMuTect_output[order(as.numeric(VCF_output$`#CHROM`),as.numeric(VCF_output$POS)), ]
  MSMuTect_output <-MSMuTect_output[!duplicated(VCF_output[, 1:5]), ]#export unique variants only!
}
VCF_output<-VCF_output[order(as.numeric(VCF_output$`#CHROM`), as.numeric(VCF_output$POS)), ]
VCF_output_unique<-VCF_output[!duplicated(VCF_output[, 1:5]), ]#export unique variants only!


setwd(outputPath)
# this function remove tunor alleles that are similar to REF (0/0)
RM_homoREF<-function(){
MSMuTect_output  <-MSMuTect_output[!VCF_output_unique$TUMOR=="0/0",]  #remove variants homozygous to REF
VCF_output_unique<-VCF_output_unique[!VCF_output_unique$TUMOR=="0/0",] #remove variants homozygous to REF
}
# this function split multiAllele to mono alleles using bcftools and save it to a new file starting with "SL_" tag.
singleLines<-function(){
    if (singleLine){
    NAME<-paste0(INPUT,'_',Sys.Date(),'_F',mode,'.vcf')
    cat(paste0(substr(Sys.time(),1,19),"  @@@ split multiAllelic loci to single lines\n"))
    system(paste0("bcftools norm -m -any ",outputPath,"/",DIRn,"/",NAME," > ",outputPath,"/",DIRn,"/","SL_",NAME))
    FNAME<-paste0(outputPath,"/",DIRn,"/","SL_",NAME)
    bcfOUT<-read.delim(FNAME,skip=9)
    bcfOUT<-bcfOUT[grepl("1",bcfOUT$TUMOR),]
    bcfOUT<-bcfOUT[!(grepl("1",bcfOUT$TUMOR)&grepl("1",bcfOUT$NORMAL)),]
    bcfOUT$TUMOR<-"0/1"
    HEADER<-readLines(FNAME,n = 10)
    write.table(HEADER,FNAME,sep = "\t",quote = F,row.names = F,col.names = F)
    write.table(bcfOUT,FNAME,sep = "\t",quote = F,row.names = F,col.names = F,append = T)
    }
}
# this function save  Files
Save_Files<-function(){
  cat(paste0("\n",substr(Sys.time(),1,19),"  @@@ Saving VCF file..."))
  write.table(VCF_header,paste0(outputPath,"/",DIRn,"/",INPUT,'_',Sys.Date(),'_F',mode,'.vcf'),sep = "\t",quote = F,row.names = F,col.names = F)
  VCF_output_unique[VCF_output_unique=="NA"]<-NA
  VCF_output_unique<-na.omit(VCF_output_unique)
  for(col in 1:ncol(VCF_output_unique)){class(VCF_output_unique[,col])<-"character"}
  write.table(VCF_output_unique,paste0(outputPath,"/",DIRn,"/",INPUT,'_',Sys.Date(),'_F',mode,'.vcf'),sep = "\t",quote = F,row.names = F,col.names = F,append = T)
  
  if(file.exists(paste0(outputPath, "/", DIRn,"/",INPUT,'_',Sys.Date(),'_F',mode,'.vcf'))) {
    cat(paste0(substr(Sys.time(), 1, 19), "  VCF file created.\n"))
  } else{
    stop("  Somthing is wrong, file not found ?!\n")
  }
  
  if (output_MSMuTect & !singleLine){
    cat(paste0(substr(Sys.time(),1,19),"  @@@ Saving filtered MSMuTect file...\n"))
    write.table(MSMuTect_output,file = paste0(outputPath,"/",DIRn,"/",gsub(pattern = "full.mut",replacement = "filtered",x = INPUT)),quote = F,col.names = T,row.names = F,sep = "\t")
    
    if(file.exists(paste0(outputPath,"/",DIRn,"/",gsub(pattern = "full.mut",replacement = "filtered",x = INPUT)))) {
      cat(paste0(substr(Sys.time(), 1, 19), "  MSMuTect.filtered file created.\n"))
    } else{
      stop("  Somthing is wrong, file not found ?!\n")
    }
  }
}

if (!OneSample) {
  cat(paste0(substr(Sys.time(),1,19),"  @@@ Normal filtering step: "))
  
  if (mode == 1) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[VCF_output_unique$NORMAL %in% c("0/0", "0|0"), ]}
    VCF_output_unique                   <-VCF_output_unique[VCF_output_unique$NORMAL %in% c("0/0", "0|0"), ]
    DIRn <- "Keep_Only_Homozygous_to_REF"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 2) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[grepl("0", VCF_output_unique$NORMAL), ]}
    VCF_output_unique <-VCF_output_unique[grepl("0", VCF_output_unique$NORMAL), ]
    DIRn <- "Filter_normals_homozygous_to_ALT"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 3) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[!grepl(pattern = "0/[1|2|3|4|5|6|7|8|9]|[1|2|3|4|5|6|7|8|9]/0", VCF_output_unique$NORMAL), ]}
    VCF_output_unique                   <-VCF_output_unique[!grepl(pattern = "0/[1|2|3|4|5|6|7|8|9]|[1|2|3|4|5|6|7|8|9]/0", VCF_output_unique$NORMAL), ]
    DIRn <- "Filter_normals_heterozygous_to_ALT"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 4) {
    DIRn <- "N_Not_Filtered"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode==5){
    RM_homoREF()
    RemoveMe <- -as.numeric(names(which(colSums(apply(VCF_output_unique[!VCF_output_unique$NORMAL %in% c("0/0","1/1") & nchar(VCF_output_unique$TUMOR)<4 ,c("NORMAL","TUMOR")],1,function(Row) strsplit(Row[2],"/")[[1]] %in% strsplit(Row[1],"/")[[1]] ))==2)))
    VCF_output_unique<-VCF_output_unique[RemoveMe,]
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[RemoveMe, ]}
    DIRn <- "Exclude_LOH"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
  }
  else if (mode==6){
    RM_homoREF()
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[sapply(strsplit(VCF_output_unique$NORMAL,"/"),function(X) X[1]==X[2]),]}
    VCF_output_unique                   <-VCF_output_unique[sapply(strsplit(VCF_output_unique$NORMAL,"/"),function(X) X[1]==X[2]),]
    DIRn <- "Exclude_normal_heterozygous"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  }  else if (mode==7){
    RM_homoREF()
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[VCF_output_unique$NORMAL %in% "0/0" & VCF_output_unique$TUMOR %in% c("0/1","1/0"), ]}
    VCF_output_unique                   <-VCF_output_unique[VCF_output_unique$NORMAL %in% "0/0" & VCF_output_unique$TUMOR %in% c("0/1","1/0"), ]
    DIRn <- "Keep_N_Homozygous_to_REF_&_T_diploid_heterozygous_tumors"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode==8){
    RM_homoREF()
    RemoveMe <- -as.numeric(names(which(colSums(apply(VCF_output_unique[!VCF_output_unique$NORMAL %in% c("0/0","1/1") & nchar(VCF_output_unique$TUMOR)<4 ,c("NORMAL","TUMOR")],1,function(Row) strsplit(Row[2],"/")[[1]] %in% strsplit(Row[1],"/")[[1]] ))==2)))
    VCF_output_unique<-VCF_output_unique[RemoveMe,]
    if (output_MSMuTect & !singleLine){MSMuTect_output <-MSMuTect_output[RemoveMe, ]}
    RemoveMe<-sapply(1:nrow(VCF_output_unique),function(ROW) all(unlist(strsplit(VCF_output_unique$TUMOR[ROW],"/")) %in% unlist(strsplit(c(VCF_output_unique$NORMAL[ROW],"0"),"/"))))   
    VCF_output_unique<-VCF_output_unique[!RemoveMe,]
    MSMuTect_output<-MSMuTect_output[!RemoveMe,]
    DIRn <- "novel_only"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
  }
  else{
    stop("\n Bad input in 'mode' argument. valid values: 1..8.\n")
  }
}

cat("\n\n\tDone.\n\n")
