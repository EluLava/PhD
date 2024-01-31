#FuniWE function
get.funiw<-function(dos){

    if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
    p<-dos@p
    het<-2*p*(1-p)
    num<-apply(as.matrix(dos),1,function(x) sum(x^2-(1+2*p)*x+2*p^2))
    den<-sum(het)
    return(list(het=den,Funi=num/den))
}

########################################################################################

#estimate FuniWE in parallel
get.funiwT<-function(bed,nloc=100000,nb.cores=30){
    
    #calculates overall Funiw for bed matrix with more than 100k SNPs
    if (dim(bed)[2]<=nloc){
        get.funiw(bed)[[2]]

    } else {
        nl<-dim(bed)[2]
        a<-split(1:nl,floor(1:nl/nloc))
        tmp_all<-parallel::mclapply(a,function(x) get.funiw(bed[,x]),mc.cores=nb.cores)
        sdens<-sum(unlist(lapply(tmp_all,function(x) x[[1]])))

        snums<-vector(length=dim(bed)[1])
        for(i in 1:length(tmp_all)){
            snums<-snums+tmp_all[[i]][[2]]*tmp_all[[i]][[1]]
        }
        snums/sdens
    }
}

########################################################################################

#defined parallelized matching function
get.matchingT <- function(bed,nloc=100000,nb.cores=4){

    #defining a new matching function (!= matching() from hierfstat)
    matching<-function(dos){
        dos <- gaston::as.matrix(dos)
        NAs<-sum(is.na(dos))
        if (NAs > 0) {
            na <- matrix(rep(1, prod(dim(dos))), ncol = ncol(dos))
            ina <- which(is.na(dos))
            na[ina] <- 0
            dos[ina] <- 1
            mNAs<-tcrossprod(na)
            Mij <- 1/2 * (1 + 1/mNAs * tcrossprod(dos - 1))
        } else {
            nl <- dim(dos)[2]
            Mij <- 1/2 * (1 + tcrossprod(dos - 1)/nl)
        }

        if (NAs==0) nl<-rep(nl,ni) else nl<-diag(mNAs)

        return(list(Mij=Mij,nl=nl))
    }

    if (class(bed)[[1]] != "bed.matrix") stop("bed must be of class bed.matrix. Exiting")

    if (dim(bed)[2]<=nloc){
        matching(bed)$Mij
    } else {
        nl<-dim(bed)[2]
        ni<-dim(bed)[1]
        a<-split(1:nl,floor(1:nl/nloc))

        tmp_all<-parallel::mclapply(a,function(x) matching(bed[,x]),mc.cores=nb.cores)

        MijT<-matrix(numeric(ni^2),ncol=ni)
        nlt<-rep(0,ni)
        for (i in 1:length(tmp_all)) {
            MijT<-MijT+tmp_all[[i]]$Mij*tmp_all[[i]]$nl
            nlt<-nlt+tmp_all[[i]]$nl
        }

        MijT/nlt

    }
}


########################################################################################


#' Interpolate genomic positions in centimorgans (CM) based on positions in base pairs (bp)
#'
#' Given a table of genomic positions in bp with their corresponding positions in CM, and a second table of SNPs
#' with positions only in bp, this function estimates the positions in CM of the SNPs in the second table based on
#' linear interpolation between the nearest SNPs in the first table.
#'
#' @param bp_table A data frame with columns 'BP' (bp position) and 'CM' (corresponding CM position)
#' @param bp_only_table A data frame with column 'BP' (bp position)
#' @return A data frame with columns 'BP' and 'CM', with the same number of rows as the input 'bp_only_table'
#' @examples
#' # Create sample data frames with some SNPs and their positions in bp and CM
#' bp_table <- data.frame(BP = c(100, 200, 300, 400, 500),
#'                        CM = c(0, 10, 20, 30, 40))
#' bp_only_table <- data.frame(BP = c(50, 150, 250, 350, 450, 550))
#' # Use the interpolate_cm function to estimate the CM positions of the SNPs in the bp_only_table
#' interp_table <- interpolate_cm(bp_table, bp_only_table)
#' # Print the result
#' print(interp_table)
#' @export
# Define a function to interpolate positions in CM based on positions in bp
interpolate_cm <- function(bp_table, bp_only_table) {

  # Sort the bp_table by bp position
  #bp_table <- bp_table[order(bp_table$BP), ]

  # Sort the bp_only_table by bp position
  #bp_only_table <- bp_only_table[order(bp_only_table$BP), ]

  # Initialize an empty vector to store interpolated CM positions
  interp_cm <- numeric(length = nrow(bp_only_table))

  # Loop through each SNP in the bp_only_table
  for (i in seq_len(nrow(bp_only_table))) {

    # Find the index of the SNP in the bp_table that comes before and after the current SNP
    before_index <- max(which(bp_table$BP < bp_only_table$BP[i]))
    after_index <- min(which(bp_table$BP >= bp_only_table$BP[i]))

    #Set value to go from physic positions to genetic distances because not the same for autosaumes and sexual chromosomes
    if("${ss}" %in% c("Super-Scaffold_13","Super-Scaffold_42")){
        val = 750000
    } else {
        val = 500000
    }

    # If the current SNP is outside the range of the bp_table, assume constant recRate for this SNP
    if (before_index==-Inf) {

      #We'll just assume position 1 (bp) is the previous one and that it's position in cM is 1/val
      DistanceBP = bp_table$BP[after_index] - 1
      DistanceCM = bp_table$CM[after_index] - 1/val

      PartInterBP = ((bp_only_table$BP[i] - 1)/DistanceBP)
      interp_cm[i] =  1/val + (DistanceCM*PartInterBP)

    #assume constant rate (with the val value) if it is after
    } else if (after_index == Inf) {

      DistanceBP = max(bp_only_table$BP) - bp_table$BP[before_index]
      DistanceCM = (DistanceBP/val)

      PartInterBP = ((bp_only_table$BP[i] - bp_table$BP[before_index])/DistanceBP)
      interp_cm[i] =  bp_table$CM[before_index] + (DistanceCM*PartInterBP)

    } else {

      DistanceBP = bp_table$BP[after_index] - bp_table$BP[before_index]
      DistanceCM = bp_table$CM[after_index] - bp_table$CM[before_index]

      PartInterBP = ((bp_only_table$BP[i] - bp_table$BP[before_index])/DistanceBP)
      interp_cm[i] =  bp_table$CM[before_index] + (DistanceCM*PartInterBP)

      ### Depreacted ###
      # Otherwise, calculate the slope and intercept of the line connecting the neighboring SNPs in the bp_table
      # slope <- (bp_table$CM[after_index] - bp_table$CM[before_index]) / (bp_table$BP[after_index] - bp_table$BP[before_index])
      # intercept <- bp_table$CM[before_index] - slope * bp_table$BP[before_index]
      # 
      # # Use the slope and intercept to interpolate the CM position of the current SNP
      # interp_cm[i] <- slope * bp_only_table$BP[i] + intercept
      ### ###
    }
  }

  # Add the interpolated CM positions to the bp_only_table and return the result
  bp_only_table$CM <- interp_cm

  #Check if the cM output is sorted
  if(is.unsorted(bp_only_table$CM)){

    print("We have a problem boss: the centiMorgans position are NOT sorted !!!")

  }

  return(bp_only_table)
}

########################################################################################

#List pairs of INDVs, used by find.unrelated function
find.names <- function(beta.mat,cutoff){
  tmp <- as.data.frame(as.matrix(which(beta.mat >= cutoff, arr.ind = T),ncol=2))
  tmp <- tmp[tmp[,1] >= tmp[,2],]
  df <- data.frame('id1'=rownames(beta.mat)[tmp[,1]],
                   'id2'=rownames(beta.mat)[tmp[,2]])
  #TEST, rm same indv pairs
  df = df[df[,1] != df[,2],]
  return(df)
}


########################################################################################

#Create new df with new pairs, used by find.unrelated function
make.ind.df <- function(pairs){
  tmp <- data.frame('id'=unique(c(pairs$id1,pairs$id2)))  # make individual df
  tmp$links <- 0 #  add link count
  for(i in 1:nrow(tmp)){
    tmp$links[i] <- sum(c(pairs$id1,pairs$id2)==tmp$id[i])
  }
  tmp <- tmp[order(tmp[,2], decreasing = T),]
  return(tmp)
}


########################################################################################

#get a list of unrelted individuals from a beta (kinship) matrix
find.unrelated <- function(b,cutoff){
  bu <- b
  # initialize pairs  
  pairs <- find.names(b,cutoff) # get pairs
  while(nrow(pairs)!=0){
    # find and remove worst candidate
    tmp <- make.ind.df(pairs)
    rm.ind <- tmp[1,1]
    bu <- bu[-match(rm.ind,rownames(bu)),-match(rm.ind,rownames(bu))]
    # recalculate pairs
    pairs <- find.names(bu,cutoff)
  }
  return(bu)
  print(timestamp())
}
