#' Reading data from Genepop
#'
#' The function \code{read.genepop} reads Genepop data files (.gen) and convert
#' them into a \linkS4class{genind} object.
#'
#' Note: \code{read.genepop} is meant for DIPLOID DATA ONLY. Haploid data with
#' the Genepop format can be read into R using \code{read.table} or
#' \code{read.csv} after removing headers and 'POP' lines, and then converted
#' using \code{\link{df2genind}}.
#'
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param ncode an integer indicating the number of characters used to code an allele.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genetix}}
#' @references Raymond M. & Rousset F, (1995). GENEPOP (version 1.2):
#' population genetics software for exact tests and ecumenicism. \emph{J.
#' Heredity}, \bold{86}:248-249 \cr
#' @keywords manip
#'
#' @export read.genepop
read.genepop <- function(file, ncode=2L, quiet=FALSE){
  ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
  if(toupper(.readExt(file)) != "GEN") stop("File extension .gen expected")
  
  if(!quiet) cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
  
  prevcall <- match.call()
  
  txt <- scan(file,sep="\n",what="character",quiet=TRUE)
  if(!quiet) cat("\nFile description: ",txt[1], "\n")
  txt <- txt[-1]
  txt <- gsub("\t", " ", txt)
  NA.char <- paste(rep("0",ncode), collapse="")
  
  ## two cases for locus names:
  ## 1) all on the same row, separated by ","
  ## 2) one per row
  ## ! spaces and tab allowed
  ## a bug was reported by S. Devillard, occuring
  ## when the two cases occur together,
  ## that is:
  ## loc1,
  ## loc2,
  ## ...
  
  
  ## new strategy (shorter): isolate the 'locus names' part and then parse it.
  locinfo.idx <- 1:(min(grep("POP",toupper(txt)))-1)
  locinfo <- txt[locinfo.idx]
  locinfo <- paste(locinfo,collapse=",")
  loc.names <- unlist(strsplit(locinfo,"([,]|[\n])+"))
  loc.names <- .rmspaces(loc.names)
  nloc <- length(loc.names)
  txt <- txt[-locinfo.idx]
  
  ## locus names have been retreived
  
  ## build the pop factor
  ## and correct the genotypes splited on more than 1 line
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
  npop <- length(pop.idx)
  ## correction for splited genotype
  ## isolated by the absence of comma on a line not containing "pop"
  nocomma <- which(! (1:length(txt)) %in% grep(",",txt))
  splited <- nocomma[which(! nocomma %in% pop.idx)]
  if(length(splited)>0){
    for(i in sort(splited,decreasing=TRUE)){
      txt[i-1] <- paste(txt[i-1],txt[i],sep=" ")
    }
    txt <- txt[-splited]
  }
  ## end correction
  
  ## reevaluate pop index
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
  
  txt[length(txt)+1] <- "POP"
  nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt)))-1
  pop <- factor(rep(1:npop,nind.bypop))
  
  txt <- txt[-c(pop.idx,length(txt))]
  
  temp <- sapply(1:length(txt),function(i) strsplit(txt[i],","))
  ## temp is a list with nind elements, first being ind. name and 2nd, genotype
  
  ind.names <- sapply(temp,function(e) e[1])
  ind.names <- .rmspaces(ind.names)
  ## individuals' name are now clean
  
  vec.genot <- sapply(temp,function(e) e[2])
  vec.genot <- .rmspaces(vec.genot)
  
  ## X is a individual x locus genotypes matrix
  X <- matrix(unlist(strsplit(vec.genot,"[[:space:]]+")),ncol=nloc,byrow=TRUE)
  
  rownames(X) <- 1:nrow(X)
  colnames(X) <- loc.names
  
  ## give right pop names
  ## beware: genepop takes the name of the last individual of a sample as this sample's name
  pop.names.idx <- cumsum(table(pop))
  pop.names <- ind.names[pop.names.idx]
  levels(pop) <- pop.names
  
  ## check that data are consistent with NCODE and ploidy=2
  if(!all(unique(nchar(X))==(ncode*2))) stop(paste("some alleles are not encoded with", ncode,
                                                   "characters\nCheck 'ncode' argument"))
  
  res <- adegenet::df2genind(X=X,pop=pop, ploidy=2, ncode=ncode, NA.char=NA.char)
  res@call <- prevcall
  
  if(!quiet) cat("\n...done.\n\n")
  
  return(res)
  
} # end read.genepop