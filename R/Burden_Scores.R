Gene_Centric_Coding_G_Star <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
                                       genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                       QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                       Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  genes <- genes_info[genes_info[,2]==chr,]

  phenotype.id <- obj_nullmodel$id_include
  ## get SNV id, position, REF, ALT (whole genome)
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant"){
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV"){
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel"){
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")

  rm(filter)
  gc()

  ### Gene
  kk <- which(genes[,1]==gene_name)

  sub_start_loc <- genes[kk,3]
  sub_end_loc <- genes[kk,4]

  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## plof
  ## Gencode_Exonic
  GENCODE.EXONIC.Category  <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

  variant.id.gene <- seqGetData(genofile, "variant.id")

  if(category == "plof"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "plof_ds"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "missense"){
    lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.missense]
  }else if(category == "disruptive_missense"){
    lof.in.ds <- ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.ds]
  }else{
    lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.synonymous]
  }

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  genotype <- Geno

  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  G <- Geno_rare
  rm(Geno_rare)
  gc()

  if(is.null(dim(G))){
    G <- matrix(G,ncol = 1)
  }

  C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)

  seqResetFilter(genofile)

  return(C)
}


Gene_Centric_Noncoding_G_Star <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                          genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                          Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){

  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  genes <- genes_info[genes_info[,2]==chr,]

  phenotype.id <- obj_nullmodel$id_include

  if(category == "downstream"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")

    rm(filter)

    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="downstream")&(SNVlist)
    variant.id.downstream <- variant.id[is.in]
    rm(GENCODE.Category)

    seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
    rm(variant.id.downstream)

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

    rm(GENCODE.Info)
    rm(variant_gene_num)

    Gene <- as.character(unlist(GENCODE.Info.split))

    rm(GENCODE.Info.split)

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else if(category == "upstream"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")

    rm(filter)

    ## upstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="upstream")&(SNVlist)
    variant.id.upstream <- variant.id[is.in]
    rm(GENCODE.Category)

    seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
    rm(variant.id.upstream)

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

    rm(GENCODE.Info)
    rm(variant_gene_num)

    Gene <- as.character(unlist(GENCODE.Info.split))
    rm(GENCODE.Info.split)

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else if(category == "UTR"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")
    rm(filter)

    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
    variant.id.UTR <- variant.id[is.in]
    rm(GENCODE.Category)

    seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
    rm(variant.id.UTR)

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

    rm(GENCODE.Info)

    # Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[seq(1,length(z),2)]))
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
    rm(GENCODE.Info.split)

    variant.id.SNV <- seqGetData(genofile, "variant.id")

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else if(category == "promoter_CAGE"){
    ## Promoter
    varid <- seqGetData(genofile, "variant.id")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

    #Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGEBvt <- CAGEAnno!=""
    CAGEidx <- which(CAGEBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEidx])
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
    ##obtain variants info
    CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
    CAGEvref <- as.character(seqGetData(genofile,"$ref"))
    CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

    rm(varid)

    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]

    dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
    dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
    dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
    dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

    seqResetFilter(genofile)
    rm(dfPromCAGEVarGene)

    ### Gene
    is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else if(category=="promoter_DHS"){
    ## Promoter
    varid <- seqGetData(genofile, "variant.id")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

    # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRsBvt <- rOCRsAnno!=""
    rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsidx])

    seqSetFilter(genofile,promGobj,intersect=TRUE)
    rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
    ## obtain variants info
    rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
    rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
    rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

    print("head(dfPromrOCRsVarGene)")
    print(head(dfPromrOCRsVarGene))

    rm(varid)

    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]

    dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
    dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
    dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
    dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

    seqResetFilter(genofile)
    rm(dfPromrOCRsVarGene)

    ### Gene
    is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)

    print("is.in")
    print(is.in)
    print("str(dfPromrOCRsVarGene.SNV)")
    print(str(dfPromrOCRsVarGene.SNV))

    variant.is.in <- variant.id.SNV[is.in]

  }else if(category=="enhancer_CAGE"){
    ## Enhancer
    varid <- seqGetData(genofile, "variant.id")

    #Now extract the GeneHancer with CAGE Signal Overlay
    genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    genehancer <- genehancerAnno!=""

    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGE <- CAGEAnno!=""
    CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
    CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

    # variants that covered by whole GeneHancer without CAGE overlap.
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

    rm(varid)

    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]

    dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
    dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
    dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
    dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

    seqResetFilter(genofile)

    rm(dfHancerCAGEVarGene)

    ### Gene
    is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else if(category=="enhancer_DHS"){
    ## Enhancer
    varid <- seqGetData(genofile, "variant.id")

    #Now extract the GeneHancer with rOCRs Signal Overlay
    genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    genehancer <- genehancerAnno!=""

    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRs <- rOCRsAnno!=""
    rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
    rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
    # variants that covered by whole GeneHancer without rOCRs overlap.

    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

    rm(varid)

    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }

    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]

    dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
    dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
    dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
    dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

    seqResetFilter(genofile)
    rm(dfHancerrOCRsVarGene)

    ### Gene
    is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

  }else{
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }

    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }

    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    variant.id <- seqGetData(genofile, "variant.id")
    rm(filter)

    ## ncRNA SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

    variant.id.ncRNA <- variant.id[is.in]
    rm(GENCODE.Category)

    seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)
    rm(variant.id.ncRNA)

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

    Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
    Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
    Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

    rm(GENCODE.Info)
    rm(GENCODE.Info.split)

    variant.id.ncRNA <- seqGetData(genofile, "variant.id")

    seqResetFilter(genofile)

    ### Gene
    is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
    is.in <- union(is.in,which(Gene_list_3==gene_name))

    variant.is.in <- variant.id.ncRNA[is.in]

  }

  print("sum(variant.is.in)")
  print(sum(variant.is.in))

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ##Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  seqResetFilter(genofile)

  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  genotype <- Geno

  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)

  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))

  print("sum(RV_label)")
  print(sum(RV_label))

  Geno_rare <- genotype$Geno[,RV_label]
  G <- Geno_rare
  rm(Geno_rare)
  gc()

  print("dim(G)")
  print(dim(G))

  if(is.null(dim(G))){
    G <- matrix(G,ncol = 1)
  }

  C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)

  return(C)
}

#' @title Extract Burden Scores from a aGDS
#' @description Extracts a non-weighted Burden Score; b = G x 1, for a given chromosome, protein-coding gene, and functional category for rare variants in either the coding region or noncoding region.
#'
#' @param region A character value specifying the region of rare variants, either Coding or Noncoding. Coding region contains 5 functional categories \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, or \code{synonymous} and noncoding region contains 8 functional categories \code{downstream}, \code{upstream}, \code{UTR}, \code{promoter_CAGE}, \code{promoter_DHS}, \code{enhancer_CAGE}, \code{enhancer_DHS}, or \code{ncRNA}.
#' @param chr A numeric value specifying chromosome in which the protein-coding gene lies in.
#' @param gene_name A character value specifying the protein-coding gene.
#' @param category A character value specifying the functional category in the protein-coding gene for which the burden score will be constructed from. When region is specified as Coding, choices include \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, or \code{synonymous} and when region is specified as Noncoding choices include \code{downstream}, \code{upstream}, \code{UTR}, \code{promoter_CAGE}, \code{promoter_DHS}, \code{enhancer_CAGE}, \code{enhancer_DHS}, or \code{ncRNA}.
#' @param genofile The opened AGDS file using SeqOpen in SeqArray.
#' @param IDs Subject IDs of individuals to construct a burden score with. Should either be a character or numeric vector matching the type from SeqGetData(genofile,"sample.id").
#' @param rare_maf_cutoff The maximum allele frequency for defining rare variants to be included in the construction of the burden score (default = 0.01).
#' @param rv_num_cutoff A numeric value specifying the cutoff of minimum number of variants to include in the construction of a burden score (default = 2).
#' @param QC_label Character value specifying the channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type Character value specifying the type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation Character value specifying the method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir Character value specifying the channel name of the annotations in the aGDS file (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog A data frame containing the name and the corresponding channel name in the aGDS file.
#' @param silent Logical: should the report of error messages be suppressed (default = FALSE).
#' @returns A column vector containing the burden score for a given chromosome, protein-coding gene, and functional category for rare variants in either the coding or noncoding region.
#' @export
Burden_Scores <- function(region = c("Coding","Noncoding"),chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                   genofile,IDs,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){

  obj_nullmodel <- list(id_include = IDs)

  if(region == "Coding"){

    if(sum(category %in% c("plof","plof_ds","missense","disruptive_missense","synonymous")) != 1){
      stop("Category must be one of plof, plof_ds, missense, disruptive_missense, or synonymous")
    }

    return(Gene_Centric_Coding_G_Star(chr = chr,gene_name = gene_name,
                                      category = category, genofile = genofile,
                                      obj_nullmodel = obj_nullmodel,rare_maf_cutoff = rare_maf_cutoff,
                                      rv_num_cutoff = rv_num_cutoff,QC_label = QC_label,
                                      variant_type = variant_type,geno_missing_imputation = geno_missing_imputation,
                                      Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, silent = silent))
  }else{
    if(sum(category %in% c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA")) != 1){
      stop("Category must be one of downstream, upstream, UTR, promoter_CAGE, promoter_DHS, enhancer_CAGE, enhancer_DHS, ncRNA")
    }
    return(Gene_Centric_Noncoding_G_Star(chr = chr,gene_name = gene_name,
                                         category = category, genofile = genofile,
                                         obj_nullmodel = obj_nullmodel,rare_maf_cutoff = rare_maf_cutoff,
                                         rv_num_cutoff = rv_num_cutoff,QC_label = QC_label,
                                         variant_type = variant_type,geno_missing_imputation = geno_missing_imputation,
                                         Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, silent = silent))
  }
}
