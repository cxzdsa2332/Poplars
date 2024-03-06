rm(list=ls())
library(topGO)
library(Biostrings)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(reshape2)
library(RColorBrewer)

load("ode_all.Rdata") # from Figure4 calculation results

check <- function(j,s,ode_all_result){
  effect_scion = mean(ode_all_result[[j]][[s]]$ode_result$Rootstock$fit$Scion)
  effect_rootstock = mean(ode_all_result[[j]][[s]]$ode_result$Scion$fit$Rootstock)
  
  name = all_names[[j]][s]  
  threthod = 0.01
  
  if (abs(effect_scion)<=threthod) {
    effect_scion = 0
  }
  
  if (abs(effect_rootstock)<=threthod) {
    effect_rootstock = 0
  }
  
  
  
  if ( effect_scion > 0 & effect_rootstock > 0){
    type = "Dove-Dove"
  } else if (effect_scion == 0 & effect_rootstock > 0) {
    type = "Zero-Dove"
  } else if (effect_scion < 0 & effect_rootstock > 0) {
    type = "Hawk-Dove"
  } else if (effect_scion > 0 & effect_rootstock == 0) {
    type = "Dove-Zero"
  } else if (effect_scion == 0 & effect_rootstock == 0) {
    type = "Zero-Zero"
  } else if (effect_scion < 0 & effect_rootstock == 0) {
    type = "Hawk-Zero"
  } else if (effect_scion > 0 & effect_rootstock < 0) {
    type = "Dove-Hawk"
  } else if (effect_scion == 0 & effect_rootstock < 0) {
    type = "Zero-Hawk"
  } else if (effect_scion < 0 & effect_rootstock < 0) {
    type = "Hawk-Hawk"
  } 
  result = list(type,name)
  return(result)
}


check_all <- function(n,ode_all_result){
  SS = data.frame(table(sapply(1:n,function(x)check(1,x,ode_all_result)[[1]])))
  ST = data.frame(table(sapply(1:n,function(x)check(2,x,ode_all_result)[[1]])))
  TS = data.frame(table(sapply(1:n,function(x)check(3,x,ode_all_result)[[1]])))
  TT = data.frame(table(sapply(1:n,function(x)check(4,x,ode_all_result)[[1]])))
  
  all_type = c("Dove-Dove","Zero-Dove","Hawk-Dove",
               "Dove-Zero","Zero-Zero","Hawk-Zero",
               "Dove-Hawk","Zero-Hawk","Hawk-Hawk")
  df_reference = data.frame(Var1 = all_type)
  df_reference = merge(df_reference,SS,by = "Var1",all.x = T)
  df_reference = merge(df_reference,ST,by = "Var1",all.x = T)
  df_reference = merge(df_reference,TS,by = "Var1",all.x = T)
  df_reference = merge(df_reference,TT,by = "Var1",all.x = T)
  colnames(df_reference) = c("type", "S/S","S/T","T/S","T/T")
  df_reference[is.na(df_reference)] = 0
  return(df_reference)
}

## note in my demonstration code, I split the whole data into 5 part and calculate each part, then merge it.
result1 = list(check_all(5000,ode_all_result = ode_all_all[[1]]),
               check_all(5000,ode_all_result = ode_all_all[[2]]),
               check_all(5000,ode_all_result = ode_all_all[[3]]),
               check_all(3000,ode_all_result = ode_all_all[[4]]),
               check_all(3252,ode_all_result = ode_all_all[[5]]))



check_gene <- function(n,list_number,ode_all_result){
  all_type = c("Dove-Dove","Zero-Dove","Hawk-Dove",
               "Dove-Zero","Zero-Zero","Hawk-Zero",
               "Dove-Hawk","Zero-Hawk","Hawk-Hawk")
  
  names = all_names[list_number]  
  
  SS = data.frame(type = sapply(1:n,function(x)check(1,x,ode_all_result)[[1]]),
                  names = names)
  colnames(SS) = c("type","name")
  #SS_all = lapply(1:9,function(x)subset(SS,SS$type==all_type[x]))
  
  
  ST = data.frame(type = sapply(1:n,function(x)check(2,x,ode_all_result)[[1]]),
                  names = names)
  #ST_all = lapply(1:9,function(x)subset(ST,ST$type==all_type[x]))
  colnames(ST) = c("type","name")
  
  TS = data.frame(type = sapply(1:n,function(x)check(3,x,ode_all_result)[[1]]),
                  names = names)
  #TS_all = lapply(1:9,function(x)subset(TS,TS$type==all_type[x]))
  colnames(TS) = c("type","name")
  
  TT = data.frame(type = sapply(1:n,function(x)check(4,x,ode_all_result)[[1]]),
                  names = names)
  #TT_all = lapply(1:9,function(x)subset(TT,TT$type==all_type[x]))
  colnames(TT) = c("type","name")
  
  result = list(SS = SS,
                ST = ST,
                TS = TS,
                TT = TT)
  
  return(result)
}


result2 = list(check_gene(5000,1,ode_all_all[[1]]),
               check_gene(5000,2,ode_all_all[[2]]),
               check_gene(5000,3,ode_all_all[[3]]),
               check_gene(3000,4,ode_all_all[[4]]),
               check_gene(3252,5,ode_all_all[[5]]))

gene_result = list(Reduce(rbind,list(result2[[1]]$SS,result2[[2]]$SS,result2[[3]]$SS,result2[[4]]$SS,result2[[5]]$SS)),
                   Reduce(rbind,list(result2[[1]]$ST,result2[[2]]$ST,result2[[3]]$ST,result2[[4]]$ST,result2[[5]]$ST)),
                   Reduce(rbind,list(result2[[1]]$TS,result2[[2]]$TS,result2[[3]]$TS,result2[[4]]$TS,result2[[5]]$TS)),
                   Reduce(rbind,list(result2[[1]]$TT,result2[[2]]$TT,result2[[3]]$TT,result2[[4]]$TT,result2[[5]]$TT)))


#SS_DOVE_DOVE = subset(gene_result[[1]],gene_result[[1]]$type=="Dove-Dove")


GOenrich <-  function(deglist,geneID2GO,useFDR = FALSE,cut=0.05, outdir = NULL,outprefix=NUL) {
  if (is.null(outdir)) {
    outdir = getwd()
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  BP.res = GOenrichsub(deglist,
                       geneID2GO,
                       'BP',
                       paste0(outdir, '/', outprefix))
  MF.res = GOenrichsub(deglist,
                       geneID2GO,
                       'MF',
                       paste0(outdir, '/', outprefix))
  CC.res = GOenrichsub(deglist,
                       geneID2GO,
                       'CC',
                       paste0(outdir, '/', outprefix))
  
  all = rbind(BP.res, MF.res, CC.res)
  all[all$Pvalue=='< 1e-30','Pvalue']=1e-31
  sigres=all[all$Pvalue < cut,]
  return(list(all,sigres))
}

GOenrichsub <- function(deglist, geneID2GO, class, output) {
  cat('Processing', class, 'class...')
  geneList <- rep(1, length(geneID2GO))
  names(geneList) <- names(geneID2GO)
  geneList[match(deglist, names(geneList))] = 0
  sampleGOdata <-
    suppressMessages(
      new(
        "topGOdata",
        nodeSize = 1,
        ontology = class,
        allGenes = geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO,
        geneSel = function(allScore) {
          return(allScore < 0.01)
        }
      )
    )
  cat('...')
  result <-
    suppressMessages(runTest(
      sampleGOdata,
      algorithm = "elim",
      statistic = "fisher"
    ))
  allRes <-
    GenTable(
      sampleGOdata,
      Pvalue = result,
      orderBy = "Pvalue",
      topNodes = attributes(result)$geneData[4],
      numChar = 1000
    )
  cat('...')
  
  gene = NULL
  GOs = sampleGOdata@graph@nodeData@data
  for (i in allRes$GO.ID) {
    e = paste0('g=names(GOs$`', i, '`$genes)')
    eval(parse(text = e))
    
    d = intersect(g, deglist)
    e = paste0('gene[\'', i, '\']=toString(d)')
    eval(parse(text = e))
  }
  allRes$gene = gene
  allRes = cbind(allRes,
                 AllDEG = length(geneList[geneList == 0]),
                 Class = class)
  cat('done!\n')
  return(allRes)
}

get_db <- function(){
  #dat = read.delim("GO.out", header = T, sep="\t")
  #GOcontent = data.frame(cbind(GID=dat$qpid,GOID=paste0('GO:', sprintf("%07d", dat$goid))))
  #geneID2GO =lapply(unique(GOcontent$GID), function (x) GOcontent$GOID[GOcontent$GID==x])
  #names(geneID2GO) = unique(GOcontent$GID)
  dat = read.delim("contig_p94_MBY_GO.txt", header = T, sep = "\t")
  dat = dat[match(rownames(df2),dat$SeqName),]
  get_list <- function(i){
    tmp = strsplit(dat[i,]$GO.IDs,"; ")
    c(unlist(lapply(tmp, substr,3,12)))
  }
  GOcontent = lapply(1:nrow(dat),function(c) get_list(c))
  #GOcontent = lapply(1:10,function(c) get_list(c))
  names(GOcontent) = unique(dat$SeqName)
  return(GOcontent)
}

geneID2GO = get_db()

#dat = read.delim(file = "contig_p94_MBY_GO.txt",sep = "\t")
#i=1;j=2
get_gene_list <- function(i,j){
  all_type = c("Dove-Dove","Zero-Dove","Hawk-Dove",
               "Dove-Zero","Zero-Zero","Hawk-Zero",
               "Dove-Hawk","Zero-Hawk","Hawk-Hawk")
  tmp = subset(gene_result[[i]],gene_result[[i]]$type==all_type[j])
  
  return(tmp$name)
}

genelist_all = lapply(1:4,function(x) lapply(1:9, function(c) get_gene_list(x,c)))

#genelist = dat$SeqName[1:1000]
#save.image(file = "GO.RData")


#GO_a <- GOenrich(genelist[[7]],geneID2GO)


GO_all <- lapply(1:4,function(x) 
  lapply(1:9,function(c) GOenrich(genelist_all[[x]][[c]],geneID2GO)))

#save.image(file = "GO.Rdata")
#dataset = GO_a
#####plot a module
dataset = GO_all[[1]][[1]]
j=1
GObubble_all <- function(j,cut = 0.05,top = 10, MainTitle = NULL) {
  require(tidyverse)
  all_type = c("Dove-Dove","Zero-Dove","Hawk-Dove",
               "Dove-Zero","Zero-Zero","Hawk-Zero",
               "Dove-Hawk","Zero-Hawk","Hawk-Hawk")
  
  conver_dataset <- function(dataset,cut = 0.05,top = 10,i){
    name = c("S/S","S/T","T/S","T/T")
    dataset$Pvalue = -log10(as.numeric(dataset$Pvalue))
    dataset = dataset %>% group_by(Class) %>% slice_max(order_by = Pvalue, n = top) %>% as.data.frame
    dataset = subset(dataset,dataset$Class == "BP")
    
    dataset$GeneRatio = dataset$Significant / dataset$AllDEG * 100
    order <- order(dataset$GeneRatio)
    dataset$Term <- factor(dataset$Term, levels = dataset$Term[order])
    dataset$type = name[i]
    dataset
  }
  dataset_all = lapply(1:4,function(c)conver_dataset(GO_all[[c]][[j]][[1]],i=c))
  dataset = Reduce(rbind,dataset_all)
  
  tmp = sapply(as.character(dataset$Term),strsplit," ")
  tmp1 = as.character(dataset$Term)[which(sapply(tmp, length)>=6)]
  tmp2 = sapply(tmp1,strsplit," ")
  a = sapply(tmp2, function(x) paste0(paste0(x[1:4],collapse = " ")," ..."))
  
  dataset$Term = as.character(dataset$Term)
  dataset$Term[which(sapply(tmp, length)>=6)] = a
  
  
  dataset_mega = lapply(1:9,function(x)lapply(1:4,function(c)conver_dataset(GO_all[[c]][[x]][[1]],i=c)))
  dataset_mega = lapply(1:9,function(x)Reduce(rbind,dataset_mega[[x]]))
  dataset_mega = Reduce(rbind,dataset_mega)
  
  xLabel = 'Count'
  yLabel = ''
  sizeLabel = 'Gene Ratio (%)'
  colorLabel = expression('-log'[10] * 'P-value')
  DotRange = c(2, 10)
  col = colorRampPalette(c("#B4FF9F","#F9F3DF", "#FFA1A1"))(100)
  col = col[15:85]
  colorbreaks = pretty(dataset_mega$Pvalue, 5)
  sizebreaks = pretty(dataset_mega$GeneRatio, 5)
  
  g <- ggplot(data = dataset, mapping = aes(Significant, Term)) +
    geom_point(aes(size = GeneRatio, color = Pvalue),
               show.legend = T,
               stroke = 1) +
    facet_grid(type ~ ., scales = "free", space = "free") +
    scale_size(
      name = sizeLabel,
      breaks = sizebreaks,
      labels = sizebreaks,
      range = DotRange,
      guide = "legend"
    ) +
    scale_color_gradientn(colours = col, breaks = colorbreaks) +
    scale_y_discrete(position = "left") +
    labs(
      title = all_type[j],
      x = xLabel,
      y = yLabel,
      size = sizeLabel,
      color = colorLabel
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(colour = "black", size = 11),
      axis.text.y = element_text(lineheight = 0.65),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = 'bold')
    )
  g
  
  return(g)
}


p_list = lapply(1:9,function(c)GObubble_all(c))

p = wrap_plots(p_list,ncol=3)
p
                
ggsave("FigS1.pdf",p, width = 30, height = 30)

