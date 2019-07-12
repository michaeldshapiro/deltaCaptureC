
## ##########################################################################
#' Downshift from DF to matrix
#' 
#' This function takes a data.frame with chr, start, end and numerical data and turns
#' it into a matrix with row names chr:start-end
#'
#' @param df This is a data frame whose first three columns are chr, start and end and
#' whose remaining columns are numerical data
#' @return A matrix of numerical data
#' @export
#' @examples
#' m = downshiftDFtoMatrix(miniSEDF)
downshiftDFtoMatrix = function(df)
{
    tags = sprintf('%s:%d-%d',df$chr,df$start,df$end)
    m = data.matrix(df[,4:ncol(df)])
    rownames(m) = tags
    return(m)
}


## ##########################################################################
#' Get the size factors for count normalization
#'
#' This function takes a data frame giving chr, start, end and count for experimental
#' replicates and returns the size factors for each of the replicates for use in
#' normalization
#'
#' @param countsDF A data frame whose first three columns are chr, start and end, and
#' whose remaining columns are count data for experimental replicates
#' @return The size factors for the columns of countsDF
#' @export
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays assay rowRanges
#' @importFrom SummarizedExperiment ranges seqnames start end
#' @importFrom SummarizedExperiment assays<- rowRanges<- mcols<-
#' @importFrom ggplot2 ggplot aes geom_col ggtitle scale_discrete_manual
#' @importFrom ggplot2 element_text theme xlab
#' @importFrom GenomicRanges GRanges seqnames isDisjoint findOverlaps
#' @importFrom GenomicRanges width intersect mcols
#' @importFrom IRanges IRanges from to width findOverlaps ranges
#' @examples
#' sf = getSizeFactorsDF(miniSEDF)
getSizeFactorsDF = function(countsDF)
{
    m = downshiftDFtoMatrix(countsDF)
    sizeFactors = estimateSizeFactorsForMatrix(m)
    return(sizeFactors)
}

## ##########################################################################
#' Get the size factors for SummarizedExperiment
#'
#' This function takes a SummarizedExperiment with an assay counts and returns
#' this object with a column sizeFactors added to its colData
#'
#' @param se A SummarizedExperiment with an assay counts
#' @return The same SummarizedExperiment with an additional column in its colData
#' giving the size factors for counts
#' @export
#' @examples
#' miniSEWithSizeFactors = getSizeFactorsSE(miniSE)
getSizeFactorsSE = function(se)
{
    colData(se)$sizeFactors = estimateSizeFactorsForMatrix(assays(se)[['counts']])
    
    return(se)
}

## ##########################################################################
#' Get normalized counts
#'
#' This function takes a SummarizedExperiment giving the the counts
#' for each replicate of the two treatments and computes and affixes
#' an assay giving the normalized version of these counts.
#'
#' @param se A SummarizedExperiment with an assay called counts giving
#'     the raw counts for each replicate of the two treatments.
#' @return A SummarizedExperiment including a an assay of the
#'     normalized counts called normalizedCounts.
#' @export
#' @examples
#' miniSENormalized = getNormalizedCountsSE(miniSE)
getNormalizedCountsSE = function(se)
{
    if(!'sizeFactors' %in% names(colData(se)))
        se = getSizeFactorsSE(se)
    
    assays(se)[['normalizedCounts']] = assays(se)[['counts']]
    for(i in seq_len(ncol(assays(se)[['normalizedCounts']])))
        assays(se)[['normalizedCounts']][,i] =
            assays(se)[['normalizedCounts']][,i] / colData(se)$sizeFactors[i]
    
    return(se)
}


## ##########################################################################
#' Make mean treatment summarized experiment:
#'
#' Get the mean normalized counts for each treatment
#'
#' This function takes a SummarizedExperiment. It looks for an assay called
#' normalizedCounts.  If this assay is missing, it creates it by normalizing
#' using the size factors. By default, it takes the mean for each value of
#' colData$treatment
#'
#' @param countsSE A SummarizedExperiment containing an assay 'counts' and
#' optionally an assay 'normalizedCounts'
#' @param byTreatment = 'treatment' This gives the column of colData to use for
#' taking averages
#' @return A SummarizedExperiment giving mean normalized counts for each value of
#' byTreatment
#' @export
#' @examples
#' meanNormalizedCountSE = getMeanNormalizedCountsSE(miniSE)
getMeanNormalizedCountsSE = function(countsSE,byTreatment='treatment')
{
    if(! 'normalizedCounts' %in% names(assays(countsSE)))
    {
        countsSE = getNormalizedCountsSE(countsSE)
    }
    
    treatments = unique(colData(countsSE)[,byTreatment])
    rho = length(treatments)
    assay = assays(countsSE)[['normalizedCounts']]
    m = matrix(0,nrow=nrow(assay),ncol=rho)
    rownames(m) = rownames(assay)
    colnames(m) = treatments
    for(tr in treatments)
    {
        idx = colData(countsSE)[,byTreatment] == tr
        ## ######################################
        ## If there's only one column, it's the answer:
        if(sum(idx) == 1)
        {
            m[,tr] = assay[,idx]
        } else { 
            m[,tr] = rowMeans(assay[,idx])
        }
    }
    colData = data.frame(treatment=treatments,
                         stringsAsFactors=FALSE)
    meanNormalizedCountsSE = SummarizedExperiment(assays=list(mean=m),
                                                  colData=colData,
                                                  rowRanges=rowRanges(countsSE))
    
    return(meanNormalizedCountsSE)
}

## ##########################################################################
#' Make delta summarized experiment:
#'
#' This function takes a SummarizedExperiment with count data and produces
#' a SummarizedExperiment of the delta track.  There should exactly two values
#' for treatment, i.e., byTreatment
#'
#' @param countsSE A summarized experiment with assay counts and optionally assay normalized
#' counts
#' @param byTreatment = 'treatment' Allows for specifying some other condition than 'treatment'
#' @return A summarized experiment with a single assay consisting of a single column, the
#' delta mean normalized counts.
#' @export
#' @examples
#' aSmallDeltaSE = getDeltaSE(miniSE)
getDeltaSE = function(countsSE,byTreatment='treatment')
{
    numTreatments = length(unique(colData(countsSE)[,byTreatment]))
    if(numTreatments != 2)
    {
        msg = paste('There should be 2 treatments. This SummarizedExperiment has',
                    numTreatments)
        stop(msg)
    }
    
    meanNormalizedCountsSE = getMeanNormalizedCountsSE(countsSE,byTreatment)
    meanCounts = assay(meanNormalizedCountsSE)
    delta = matrix(meanCounts[,1] - meanCounts[,2],ncol=1)
    colData = data.frame(delta=sprintf('%s - %s',
                                       as.character(colData(meanNormalizedCountsSE)$treatment[1]),
                                       as.character(colData(meanNormalizedCountsSE)$treatment[2])),
                         stringsAsFactors=FALSE)
    deltaSE = SummarizedExperiment(assay=list(delta=delta),
                                   colData=colData)
    
    rowRanges(deltaSE) = rowRanges(meanNormalizedCountsSE)
    
    return(deltaSE)
}

## ##########################################################################
#' Get the binning factors for one set of GRanges into another
#'
#' This function takes two GRanges, one representing a set of bins and the other representing
#' data to be pro-rated over those bins and returns a data frame giving the overlaps, various widths
#' and the fractions for pro-rating scores
#'
#' @param bins a set of GRanges to be used for binning data.  
#' @param gr the GRanges of the data to be binned
#' @param checkDisjoint = FALSE if this is TRUE it will check to see that the ranges in each of
#' bins and gr are disjoint
#' @return A data frame giving index pairs for the intersections, widths of the intersections
#' and the fraction of each gr range meeting each bin
#' @export
#' @examples
#' overlapWeights = getOverlapWeights(weightsExampleBins,weightsExampleGr)
getOverlapWeights = function(bins,gr,checkDisjoint=FALSE)
{
    if(checkDisjoint)
    {
        if(! isDisjoint(bins))
            stop('Attempting to bin data into overlapping bins.')
    }
    
    meets = findOverlaps(bins,gr)
    
    L = length(meets)
    meetsDF = data.frame(from=numeric(L),
                         to=numeric(L),
                         toWidth=numeric(L),
                         overlap=numeric(L),
                         fraction=numeric(L),
                         stringsAsFactors=FALSE)
    
    meetsDF$from = from(meets)
    meetsDF$to = to(meets)
    meetsDF$fromWidth = width(bins[from(meets)])
    meetsDF$toWidth = width(gr[to(meets)])
    
    for(i in seq_len(L))
        meetsDF$overlap[i] = width(intersect(bins[meetsDF$from[i]],gr[meetsDF$to[i]]))
    
    
    meetsDF$fraction = meetsDF$overlap / meetsDF$toWidth
    meetsDF = meetsDF[,c('from','to','overlap','fraction')]
    
    return(meetsDF)                   
}



## ##########################################################################
#' Bin a Summarized experiment into a set of bins given by a GRanges object
#'
#' This function takes a set of bins given by a GRanges object and a
#' RangedSummarizedExperiment and produces a new RangedSummarizedExperiment
#' with the bins as its rowRanges
#'
#' @param bins a GRanges object whose ranges should be disjoint
#' @param se a RangedSummarizedExperiment
#' @param checkDisjoint = FALSE if set to true will check that the bins
#'     are disjoint
#' @return a RangedSummarizedExperiment
#' @export
#' @examples
#' binnedSummarizedExperiment = binSummarizedExperiment(smallSetOfSmallBins,smallerDeltaSE)
binSummarizedExperiment = function(bins,se,checkDisjoint=FALSE)
{
    if(checkDisjoint)
    {
        if(! isDisjoint(bins))
            stop('Attempting to bin data into overlapping bins.') 
    }
    overlapsDF = getOverlapWeights(bins,rowRanges(se))
    
    binnedAssays = list()
    
    rows = length(bins)
    cols = ncol(assays(se)[[1]])
    
    for(n in names(assays(se)))
    {
        m = matrix(0,nrow=rows,ncol=cols)
        val = nrow(overlapsDF)
        for(i in seq_len(val))
        {
            whichRow = overlapsDF$from[i]
            m[whichRow,] = m[whichRow,] +
                overlapsDF$fraction[i] * assays(se)[[n]][overlapsDF$to[i],]
        }
        binnedAssays[[n]] = m
    }
    binnedSE = SummarizedExperiment(assays=binnedAssays,
                                    colData=colData(se),
                                    rowRanges=bins)
    return(binnedSE)
}

## ##########################################################################
#' Rebin a SummarizedExperiment to a multiple of its bin width
#'
#' This is a faster way of rebinning when the old bins are consecutive
#' and constant width and the new bins are to be a multiple of that width
#'
#' @param se a RangedSummarizedExperiment to be rebinned
#' @param multiple the factor by which to fatten the bins
#' @param deleteShort = FALSE when set to true if the final bin is short
#'     it will be deleted
#' @return a RangedSummarizedExperiment
#' @export
#' @examples
#' rebinnedSummarizedExperiment = rebinToMultiple(binnedSummarizedExperiment,10)
rebinToMultiple = function(se,multiple,deleteShort=FALSE)
{
    ## ######################################
    ## Bins should be constant width:
    W = width(ranges(se))
    if(length(unique(W)) != 1)
        stop('Bins should be of equal length.')
    W = W[1]
    
    ## ######################################
    ## Bins should be on one chromosome:
    if(length(unique(seqnames(se))) != 1)
        stop('Data should be binned into bins on a single chromosome.')
    
    ## ######################################
    ## Bins should be consecutive:
    n = length(se)
    theJumps = start(se)[2:n] - end(se)[seq_len(n-1)]
    theJumps = unique(theJumps)
    if(length(theJumps) !=1 |
       theJumps != 1)
        stop('Data should be binned into bins which are consecutive.')

    
    if(deleteShort)
    {
        n = multiple * floor(length(se) / multiple)
        se = se[seq_len(n)]
    }
    
    startRangeIdx = seq(1,length(se),by=multiple)
    endRangeIdx = startRangeIdx + multiple - 1
    N = length(startRangeIdx)
    endRangeIdx[N] = min(endRangeIdx[N],length(se))
    
    gr = GRanges(seqnames=seqnames(se)[1],
                 IRanges(start=start(se)[startRangeIdx],
                         end=end(se)[endRangeIdx]))
    
    binnedAssays = list()
    rows = N
    cols = ncol(assay(se))
    for(n in names(assays(se)))
    {
        m = matrix(0,nrow=rows,ncol=cols)
        for(i in seq_len(N))
            m[i,] = sum(assays(se)[[n]][(startRangeIdx[i]:endRangeIdx[i]),])
        
        binnedAssays[[n]] = m
    }
    binnedSE = SummarizedExperiment(assays=binnedAssays,
                                    colData=colData(se),
                                    rowRanges=gr)
    
    return(binnedSE)
}


## ##########################################################################
#' Generate permutation for permutation testing
#'
#' This function takes a set of row ranges and an inner region and
#' generates a permutation which is symmetric on the inner region and
#' arbitrary on the remainder
#'
#' @param gr a GRanges object which should be ordered
#' @param innerRegion a GRanges object which should be a single interval
#' @return a permutation of 1:length(gr)
#' @export
#' @examples
#' permutations = generatePermutation(smallBins,viewpointRegion)
generatePermutation = function(gr,innerRegion)
{
    ## ######################################    
    meets = findOverlaps(gr,innerRegion)
    meetsIdx = meets@from
    
    bigIdx = seq_len(length(gr))
    m = min(bigIdx[meetsIdx])
    n = max(bigIdx[meetsIdx])
    permutation = bigIdx
    ## ######################################
    ## Permute the inner region:
    val = ceiling(n-m) / 2
    for(i in seq_len(val))
    {
        if(stats::runif(1) < .5)
        {
            I = m + i
            J = n - i
            permutation[I] = J
            permutation[J] = I
        }
    }
    ## ######################################
    ## Permute the remainder:
    a = seq_len(m-1)
    b = (n+1):length(permutation)
    outer = c(a,b)
    outer = sample(outer,length(outer),replace=FALSE)
    permutation[a] = outer[a]
    permutation[b] = outer[m:length(outer)]
    
    return(permutation)
}

## ##########################################################################
#' Get the runs and their values
#'
#' This function finds the runs of consecutive ranges in which the
#' sign of the data does not change.  It returns a GRanges object
#' containing the contiguous ranges and the weighted sum of  data in
#' each.
#'
#' @param se a SummarizedExperiment whose first assay has a column named
#'     colName. Typically this will be a one-column matrix with delta.
#' @param innerRegion a Granges object defining the region surrounding
#'     the viewpoint to be excluded from run total calculations
#' @param colName defaults to 'delta'
#'
#' @return a GRanges object giving the contiguous regions and their
#'     respective sums
#' @export
#' @examples
#' runTotals = getRunTotals(binnedDeltaSE,viewpointRegion)
getRunTotals = function(se,innerRegion,colName='delta')
{
    ## ######################################
    ## We turn the SummarizedExperiment into a GRanges object:
    gr = rowRanges(se)
    mcols(gr)[,colName] = assay(se)[,1]
    
    
    
    ## ######################################
    ## Subset to miss the inner region:
    idx = end(gr) < min(start(innerRegion))
    firstHalf = gr[idx]
    idx = max(end(innerRegion)) < start(gr)
    secondHalf = gr[idx]
    

    
    return(c(.getRunsAndTotals(firstHalf,colName),
             .getRunsAndTotals(secondHalf,colName)))
}

## ##########################################################################
#' A helper function for getRunTotals
#'
#' This takes a GRanges object for binneed data and a column name
#' designating where to find the relevant data in the mcols and
#' returns a GRanges giving the consecutive runs of constant sign and
#' their run totals.  It is not exported.
#'
#' @param gr a GRanges object whose mcols gives the relevant binned
#' data
#' @param colName This designates the column in mcols with the relevant
#' data
#' @return a GRanges object giving the contiguous regions and their
#' respective sums.
.getRunsAndTotals = function(gr,colName)
{
    ## ######################################
    ## Pick up the first value, first sign, first total:
    runTotals = c()
    start = c()
    end = c()
    finger = 1
    values = mcols(gr)[,colName]
    currentValue = values[finger]
    start = c(start,start(gr)[finger])
    currentTotal = currentValue
    if(currentValue >= 0)
        currentSign = 1
    else
        currentSign = -1
    
    ## ######################################
    ## Iterate through.  At each step we either update
    ## the total or save the total, change the sign, and
    ## start a new total:
    while(finger < length(gr))
    {
        
        finger = finger + 1
        currentValue = values[finger]
        if(currentValue * currentSign >= 0)
        {
            currentTotal = currentTotal + currentValue
        }
        else
        {
            runTotals = c(runTotals,currentTotal)
            end = c(end,end(gr)[finger-1])
            start = c(start,start(gr)[finger])
            currentSign = - currentSign
            currentValue = values[finger]
            currentTotal = currentValue
        }
    }
    ## ######################################
    ## Pick up the last running total:
    runTotals = c(runTotals,currentTotal)
    end = c(end,end(gr)[finger])
    
    chr = rep(unique(seqnames(gr)),length(start))
    runTotals = GRanges(seqnames=chr,
                        IRanges(start,end),
                        runTotals=runTotals)
    
    return(runTotals)
    
}

## ##########################################################################
#' Get the lopsidedness statistic
#'
#' This function looks at the sidedness around the viewpoint and
#'  returns the absolute value of the difference between the sum of
#'  the values before and after the viewpoint inside the viewpoint
#'  region.
#'
#' @param se a SummerizedExperiment giving the delta or permuted delta
#' @param viewpointRegion the region around the viewpoint in which to
#' investigate lopsidedness
#' @param colName defaults to 'delta'
#' @return the lopsidedness around the viewpointMid in the
#'     viewpointRegion
#' @export
#' @examples
#' lopsidedness = getLopsidedness(binnedDeltaSE,viewpointRegion)
getLopsidedness = function(se,viewpointRegion,colName='delta')
{
    ## ######################################
    ## We'll want the midpoint of the viewpointRegion:
    viewpointMid = floor((start(viewpointRegion) +
                          end(viewpointRegion)) / 2)
    
    ## ######################################
    ## We turn the SummarizedExperiment into a GRanges object:
    gr = rowRanges(se)
    mcols(gr)[,colName] = assay(se)[,1]
    
    ## ######################################
    ## We want to restrict to the viewpointRegion:
    window = findOverlaps(viewpointRegion,gr)
    windowBins = to(window)
    gr = gr[windowBins]
    
    ## ######################################
    ## We want to split this below and above the midpoint:
    belowIdx = end(gr) < viewpointMid
    aboveIdx = viewpointMid < start(gr)
    below = gr[belowIdx]
    above = gr[aboveIdx]
    
    ## ######################################
    ## We want the same number of bins below and above and
    ## if we need to discard bins, we want to do this furthest
    ## from the midpoint:
    b = length(below)
    a = length(above)
    use = min(b,a)
    
    below = below[(b-use+1):b]
    above = above[seq_len(use)]
    
    lopsidedness = abs(sum(mcols(below)[,colName]) -
                       sum(mcols(above)[,colName]))
    
    return(lopsidedness)
}

## ##########################################################################
#' Get the distribution of run and lopsidedness statistics
#'
#' @param scrambledDeltas a list of rebinned (i.e., to large bin size)
#'     of scrambled deltas
#' @param viewpointRegion a GRanges object giving the region that is
#'     reserved for lopsidedness
#' @param colName = 'delta' 
#' @return a Nx4 matrix giving the min, max, max(abs(min),abs(max))
#'     and lopsidedness for the run totals in the list of scrambled
#'     deltas. 
#' @export
getRunAndLopsidednessStatistics =
    function(scrambledDeltas,viewpointRegion,colName='delta')
{
    ## ######################################
    ## We'll want the midpoint of the viewpointRegion:
    viewpointMid = floor((start(viewpointRegion) +
                          end(viewpointRegion)) / 2)
    
    scrambledRuns = lapply(scrambledDeltas,function(x)
        return(getRunTotals(x,viewpointRegion)))
    m = getRunStatisticsDist(scrambledRuns)
    lopsidedness = unlist(lapply(scrambledDeltas,function(x)
        return(getLopsidedness(x,viewpointRegion,colName='delta'))))
    m = cbind(m,lopsidedness)
    colnames(m)[4] = 'lopsidedness'
    
    return(m)
}

## ##########################################################################
#' This takes a list of (scrambled) runs and returns their run statistics
#'
#' This function takes a list of (scrambled) runs and extracts their
#' run totals as a matrix with colnames 'min','max' and 'abs', the
#' latter being the max of the absolute values of the previous two
#' 
#' @param runTotalsList this is a list whose members are GRanges
#'     objects giving the consecutive runs and their totals
#' @return a Nx3 matrix giving the min, max and max(abs(min),abs(max))
#'     run totals
#' @export
getRunStatisticsDist = function(runTotalsList)
{
    m = matrix(0,ncol=3,nrow=length(runTotalsList))
    colnames(m) = c('min','max','abs')
    val = length(runTotalsList)
    for(i in seq_len(val))
        m[i,] = getRunStatistics(runTotalsList[[i]])
    
    return(m)
}

## ##########################################################################
#' This function is called by getRunsStatisticsDist on the individual elements
#' of a list of scrambled runs.
#'
#' This is a helper function.  Currently not exported.
#'
#' @param runTotals is a GRanges object giving the consecutive runs
#'     and their totals. 
#' @return a vector of the min, max and absolute value of the min and
#'     max for the run totals. 
getRunStatistics = function(runTotals)
{
    totals = mcols(runTotals)[,'runTotals']
    return(c('min'=min(totals),
             'max'=max(totals),
             'abs'=max(abs(totals))))
}

## ##########################################################################
#' This function returns the significance levels for min, max, "abs"
#' and lopsidedness.
#'
#' Given an Nx4 matrix with columns 'min','max','abs' and
#' 'lopsidededness', this function returns the cutoff levels for a given
#' pValue.
#'
#' @param runStats a matrix with columns 'min','max','abs' and
#' 'lopsidededness'
#' @param p =.05 the desired p-value
#' @return a vector with cutoff values
#' @export
#' @examples
#' dimnames = list(c(),c('min','max','abs','lopsidedness'))
#' m = 10 * (matrix(runif(400),ncol=4,dimnames=dimnames) - 0.5)
#' cutoffs = getPValueCutoff(m,.05)
getPValueCutoff = function(runStats,p=.05)
{
    N = nrow(runStats)
    cutAt = ceiling((1-p) * N)
    mins = runStats[,'min']
    mins = mins[order(-mins)]
    minCutoff = mins[cutAt]
    
    maxs = runStats[,'max']
    maxs = maxs[order(maxs)]
    maxCutoff = maxs[cutAt]
    
    abss = runStats[,'abs']
    abss = abss[order(abss)]
    absCutoff = abss[cutAt]
    
    lops = runStats[,'lopsidedness']
    lops = lops[order(lops)]
    lopsCutoff = lops[cutAt]
    
    return(c('min'=minCutoff,
             'max'=maxCutoff,
             'abs'=absCutoff,
             'lopsidedness'=lopsCutoff))
}


## ##########################################################################
#' Get ths significant regions from delta data
#'
#' This function takes delta data as a SummarizedExperiment and
#' required ancillary data and returns a GenomicRanges object whose
#' mcols indicate the significant regions.
#'
#'
#' @param deltaSE a ranged summarized experiment with a one-column
#'     assay giving the delta mean count
#' @param regionOfInterest a GenomicRanges object specifying the
#'     region  of interest
#' @param viewpointRegion the region withheld from arbitrary
#'     permutation
#' @param smallBinSize size to bin original data to for permutation
#' @param bigBinSize size to bin data to for significance testing.  Must be a multiple of smallBinSize
#' @param numPermutations = 1000 the number of permutations to be used for
#'     permutation testing
#' @param pValue the desired significance level
#' @return a GRanges object giving the bigBin binning of region of
#'     interest whose mcols gives the values of delta and logicals
#'     telling whether the bin is in the viewpoint regsion and whether
#'     it rises to statistical significance
#' @export
getSignificantRegions = function(deltaSE,
                                 regionOfInterest,
                                 viewpointRegion,
                                 smallBinSize,
                                 bigBinSize,
                                 numPermutations=1000,
                                 pValue=0.05)
{
    ## ##########################################################################
    ## Trim to region of interest:
    chr = as.character(seqnames(regionOfInterest)[1])
    fromHere = min(start(regionOfInterest))
    toHere = max(end(regionOfInterest))
    idx = (as.character(seqnames(deltaSE)) == chr &
           fromHere <= start(deltaSE) &
           end(deltaSE) <= toHere)
    deltaSE = deltaSE[idx,]
    
    ## ##########################################################################
    ## Bin to small bin size.  These give the start and end of each bin:
    start = seq(from=fromHere,to=toHere,by=smallBinSize)
    end = start - 1 + smallBinSize
    gr = GRanges(seqnames=chr,
                 IRanges(start,end))
    
    
    binnedDeltaSE = binSummarizedExperiment(gr,deltaSE)
    
    ## ##########################################################################
    ## Bin to big bin size:
    lambda = bigBinSize / smallBinSize
    bigBinDeltaSE = rebinToMultiple(binnedDeltaSE,lambda)
    
    ## ##########################################################################
    ## Get real runs:
    realRunTotals = getRunTotals(bigBinDeltaSE,viewpointRegion)
    mid = floor((start(viewpointRegion) +
                 end(viewpointRegion)) / 2)
    realLopsidedness = getLopsidedness(bigBinDeltaSE,viewpointRegion,mid)
    
    scrambledDeltas = list()
    scrambledRuns = list()
    rowRanges=rowRanges(binnedDeltaSE)
    colData=colData(binnedDeltaSE)
    for(j in seq_len(numPermutations))
    {
        permutation = generatePermutation(binnedDeltaSE,viewpointRegion)
        m = matrix(assay(binnedDeltaSE)[permutation,],ncol=1)
        
        scrambledDeltaSE = SummarizedExperiment(colData=colData,
                                                assays=list('delta'=m),
                                                rowRanges=rowRanges)
        ## ##########################################################################
        scrambledDeltas[[j]] = rebinToMultiple(scrambledDeltaSE,lambda)
        scrambledRuns[[j]] = getRunTotals(scrambledDeltaSE,viewpointRegion)
    }
    
    ## ##########################################################################
    ## Get run stats:
    stats = getRunAndLopsidednessStatistics(scrambledDeltas,viewpointRegion,mid)
    runStats = getRunStatisticsDist(scrambledRuns)
    
    ## ##########################################################################
    ## Determine significant levels:
    cutoffs = getPValueCutoff(stats,pValue)
    
    ## ##########################################################################
    ## Extract the results:
    answers = rowRanges(bigBinDeltaSE)
    
    mcols = data.frame(delta=assay(bigBinDeltaSE),
                       isViewpoint=FALSE,
                       lopsidednessIsSignificant=FALSE,
                       isSignificantByMin=FALSE,
                       isSignificantByMax=FALSE,
                       isSignificantByAbs=FALSE,
                       stringsAsFactors=FALSE)
    
    meets = findOverlaps(answers,viewpointRegion)
    mcols$isViewpoint[from(meets)] = TRUE
    
    ## ##########################################################################
    ## Is lopsidedness significant?
    if(realLopsidedness > cutoffs['lopsidedness'])
        mcols$lopsidednessIsSignificant[mcols$isViewpoint] = TRUE
    
    ## ##########################################################################
    ## What about runs?
    val = length(realRunTotals)
    for(i in seq_len(val))
    {
        if(realRunTotals$runTotals[i] < cutoffs['min'])
        {
            meets = findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByMin[from(meets)] = TRUE
        }
        
        if(realRunTotals$runTotals[i] > cutoffs['max'])
        {
            meets = findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByMax[from(meets)] = TRUE
        }
        
        
        if(abs(realRunTotals$runTotals[i]) > cutoffs['abs'])
        {
            meets = findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByAbs[from(meets)] = TRUE
        }
    }
    mcols(answers) = mcols
    
    return(answers)
}


## ##########################################################################
## ##########################################################################
#' This produces a plot of the region of interest showing regions of significance.
#'
#' This function takes a input the GRanges object produced by
#' getSignificant regions and produces a ggplot of significant features
#'
#'
#' @param significantRegions a GRanges object as produced by
#'     getSignificantRegions
#' @param significanceType = 'abs' a variable indicating whether to plot
#'     significance according to min, max or abs.
#' @param title a title for the plot
#' @param xLabel = 'viewpoint' supplies an xlabel
#' @param legend = TRUE whether or not to show the legend
#' @return a ggplot object
#' @export
#' @examples
#' plotOfSignificantRegions = plotSignificantRegions(significantRegions)
    plotSignificantRegions =
    function(significantRegions,significanceType='abs',title='Significant Regions',xLabel='viewpoint',legend=TRUE)
{
    sigNames = c(min="isSignificantByMin",
                 max="isSignificantByMax",
                 abs="isSignificantByAbs")
    
    df = data.frame(significantRegions,
                    stringsAsFactors=FALSE)
    
    df$significance = 'not significant'
    df$significance[df$lopsidednessIsSignificant] =
        'lopsidedness is significant'
    df$significance[df[,sigNames[significanceType]]] = 'significant'
    
    df$ctr = (df$start + df$end) / 2
    
    breaks = c('not significant','lopsidedness is significant','significant')
    values = c('#F97465','#01BB34','#6998FF')
    names(values) = breaks
    
    g = ggplot(df,aes(x=ctr,y=delta,fill=significance)) +
        geom_col() +
        ggtitle(title) +
        scale_discrete_manual(breaks=breaks,
                              values=values,
                              aesthetics=c('color','fill')) +
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=16,face="bold")) +
        xlab(xLabel)
    
    if(! legend)
        g = g + theme(legend.position='none')
    
    return(g)
}

