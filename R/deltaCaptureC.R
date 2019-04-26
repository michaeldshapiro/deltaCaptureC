
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
#' @examples
#' sf = getSizeFactorsDF(miniSEDF)
getSizeFactorsDF = function(countsDF)
{
    m = downshiftDFtoMatrix(countsDF)
    sizeFactors = DESeq2::estimateSizeFactorsForMatrix(m)
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
    SummarizedExperiment::colData(se)$sizeFactors = DESeq2::estimateSizeFactorsForMatrix(SummarizedExperiment::assays(se)[['counts']])
    
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
    if(!'sizeFactors' %in% names(SummarizedExperiment::colData(se)))
        se = getSizeFactorsSE(se)
    
    SummarizedExperiment::assays(se)[['normalizedCounts']] = SummarizedExperiment::assays(se)[['counts']]
    for(i in seq_len(ncol(SummarizedExperiment::assays(se)[['normalizedCounts']])))
        SummarizedExperiment::assays(se)[['normalizedCounts']][,i] =
            SummarizedExperiment::assays(se)[['normalizedCounts']][,i] / SummarizedExperiment::colData(se)$sizeFactors[i]
    
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
    if(! 'normalizedCounts' %in% names(SummarizedExperiment::assays(countsSE)))
    {
        countsSE = getNormalizedCountsSE(countsSE)
    }
    
    treatments = unique(SummarizedExperiment::colData(countsSE)[,byTreatment])
    rho = length(treatments)
    assay = SummarizedExperiment::assays(countsSE)[['normalizedCounts']]
    m = matrix(0,nrow=nrow(assay),ncol=rho)
    rownames(m) = rownames(assay)
    colnames(m) = treatments
    for(tr in treatments)
    {
        idx = SummarizedExperiment::colData(countsSE)[,byTreatment] == tr
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
    meanNormalizedCountsSE = SummarizedExperiment::SummarizedExperiment(assays=list(mean=m),
                                                                        colData=colData,
                                                                        rowRanges=SummarizedExperiment::rowRanges(countsSE))
    
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
    stopifnot(length(unique(SummarizedExperiment::colData(countsSE)[,byTreatment])) == 2)
    
    meanNormalizedCountsSE = getMeanNormalizedCountsSE(countsSE,byTreatment)
    meanCounts = SummarizedExperiment::assay(meanNormalizedCountsSE)
    delta = matrix(meanCounts[,1] - meanCounts[,2],ncol=1)
    colData = data.frame(delta=sprintf('%s - %s',
                                       as.character(SummarizedExperiment::colData(meanNormalizedCountsSE)$treatment[1]),
                                       as.character(SummarizedExperiment::colData(meanNormalizedCountsSE)$treatment[2])),
                         stringsAsFactors=FALSE)
    deltaSE = SummarizedExperiment::SummarizedExperiment(assay=list(delta=delta),
                                                         colData=colData)
    
    SummarizedExperiment::rowRanges(deltaSE) = SummarizedExperiment::rowRanges(meanNormalizedCountsSE)
    
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
#' overlapWeights = getOverlapWeights(smallBins[1:20],GenomicRanges::GRanges(GenomicRanges::seqnames(miniDeltaSE),GenomicRanges::ranges(miniDeltaSE))[1:40])
getOverlapWeights = function(bins,gr,checkDisjoint=FALSE)
{
    if(checkDisjoint)
    {
        stopifnot(GenomicRanges::isDisjoint(bins))
    }
    
    meets = IRanges::findOverlaps(bins,gr)
    
    L = length(meets)
    meetsDF = data.frame(from=numeric(L),
                         to=numeric(L),
                         toWidth=numeric(L),
                         overlap=numeric(L),
                         fraction=numeric(L),
                         stringsAsFactors=FALSE)
    
    meetsDF$from = IRanges::from(meets)
    meetsDF$to = IRanges::to(meets)
    meetsDF$fromWidth = GenomicRanges::width(bins[IRanges::from(meets)])
    meetsDF$toWidth = GenomicRanges::width(gr[IRanges::to(meets)])

    for(i in seq_len(L))
        meetsDF$overlap[i] = IRanges::width(GenomicRanges::intersect(bins[meetsDF$from[i]],gr[meetsDF$to[i]]))

    
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
        stopifnot(GenomicRanges::isDisjoint(bins))
    }
    overlapsDF = getOverlapWeights(bins,SummarizedExperiment::rowRanges(se))
    
    binnedAssays = list()
    
    rows = length(bins)
    cols = ncol(SummarizedExperiment::assays(se)[[1]])
    
    for(n in names(SummarizedExperiment::assays(se)))
    {
        m = matrix(0,nrow=rows,ncol=cols)
        for(i in seq_len(nrow(overlapsDF)))
        {
            whichRow = overlapsDF$from[i]
            m[whichRow,] = m[whichRow,] +
                overlapsDF$fraction[i] * SummarizedExperiment::assays(se)[[n]][overlapsDF$to[i],]
        }
        binnedAssays[[n]] = m
    }
    binnedSE = SummarizedExperiment::SummarizedExperiment(assays=binnedAssays,
                                                          colData=SummarizedExperiment::colData(se),
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
    W = IRanges::width(SummarizedExperiment::ranges(se))
    stopifnot(length(unique(W)) == 1)
    W = W[1]
    
    ## ######################################
    ## Bins should be on one chromosome:
    stopifnot(length(unique(SummarizedExperiment::seqnames(se))) == 1)
    
    ## ######################################
    ## Bins should be consecutive:
    n = length(se)
    theJumps = SummarizedExperiment::start(se)[2:n] - SummarizedExperiment::end(se)[seq_len(n-1)]
    theJumps = unique(theJumps)
    stopifnot(length(theJumps) == 1 &
              theJumps == 1)
    
    if(deleteShort)
    {
        n = multiple * floor(length(se) / multiple)
        se = se[seq_len(n)]
    }
    
    startRangeIdx = seq(1,length(se),by=multiple)
    endRangeIdx = startRangeIdx + multiple - 1
    N = length(startRangeIdx)
    endRangeIdx[N] = min(endRangeIdx[N],length(se))
    
    gr = GenomicRanges::GRanges(seqnames=SummarizedExperiment::seqnames(se)[1],
                                IRanges::IRanges(start=SummarizedExperiment::start(se)[startRangeIdx],
                                                 end=SummarizedExperiment::end(se)[endRangeIdx]))
    
    binnedAssays = list()
    rows = N
    cols = ncol(SummarizedExperiment::assay(se))
    for(n in names(SummarizedExperiment::assays(se)))
    {
        m = matrix(0,nrow=rows,ncol=cols)
        for(i in seq_len(N))
            m[i,] = sum(SummarizedExperiment::assays(se)[[n]][(startRangeIdx[i]:endRangeIdx[i]),])
        
        binnedAssays[[n]] = m
    }
    binnedSE = SummarizedExperiment::SummarizedExperiment(assays=binnedAssays,
                                                          colData=SummarizedExperiment::colData(se),
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
    meets = IRanges::findOverlaps(gr,innerRegion)
    meetsIdx = meets@from
    
    bigIdx = seq_len(length(gr))
    m = min(bigIdx[meetsIdx])
    n = max(bigIdx[meetsIdx])
    permutation = bigIdx
    ## ######################################
    ## Permute the inner region:
    for(i in seq_len(ceiling((n-m)/2)))
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
#' @return a GRanges object giving the contiguous region and their
#'     respective sums
#' @export
#' @examples
#' runTotals = getRunTotals(binnedDeltaSE,viewpointRegion)
getRunTotals = function(se,innerRegion,colName='delta')
{
    ## ######################################
    ## We turn the SummarizedExperiment into a GRanges object:
    gr = SummarizedExperiment::rowRanges(se)
    GenomicRanges::mcols(gr)[,colName] = SummarizedExperiment::assay(se)[,1]
    
    
    
    ## ######################################
    ## Subset to miss the inner region:
    idx = GenomicRanges::end(gr) < min(GenomicRanges::start(innerRegion))
    firstHalf = gr[idx]
    idx = max(GenomicRanges::end(innerRegion)) < GenomicRanges::start(gr)
    secondHalf = gr[idx]
    
    getRunsAndTotals = function(gr,colName)
    {
        ## ######################################
        ## Pick up the first value, first sign, first total:
        runTotals = c()
        start = c()
        end = c()
        finger = 1
        values = GenomicRanges::mcols(gr)[,colName]
        currentValue = values[finger]
        start = c(start,GenomicRanges::start(gr)[finger])
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
                end = c(end,GenomicRanges::end(gr)[finger-1])
                start = c(start,GenomicRanges::start(gr)[finger])
                currentSign = - currentSign
                currentValue = values[finger]
                currentTotal = currentValue
            }
        }
        ## ######################################
        ## Pick up the last running total:
        runTotals = c(runTotals,currentTotal)
        end = c(end,GenomicRanges::end(gr)[finger])
        
        chr = rep(unique(GenomicRanges::seqnames(gr)),length(start))
        runTotals = GenomicRanges::GRanges(seqnames=chr,
                                           IRanges::IRanges(start,end),
                                           runTotals=runTotals)
        
        return(runTotals)
        
    }
    return(c(getRunsAndTotals(firstHalf,colName),
             getRunsAndTotals(secondHalf,colName)))
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
    viewpointMid = floor((GenomicRanges::start(viewpointRegion) +
                          GenomicRanges::end(viewpointRegion)) / 2)
    
    ## ######################################
    ## We turn the SummarizedExperiment into a GRanges object:
    gr = SummarizedExperiment::rowRanges(se)
    GenomicRanges::mcols(gr)[,colName] = SummarizedExperiment::assay(se)[,1]
    
    ## ######################################
    ## We want to restrict to the viewpointRegion:
    window = IRanges::findOverlaps(viewpointRegion,gr)
    windowBins = IRanges::to(window)
    gr = gr[windowBins]
    
    ## ######################################
    ## We want to split this below and above the midpoint:
    belowIdx = GenomicRanges::end(gr) < viewpointMid
    aboveIdx = viewpointMid < GenomicRanges::start(gr)
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
    
    lopsidedness = abs(sum(GenomicRanges::mcols(below)[,colName]) -
                       sum(GenomicRanges::mcols(above)[,colName]))
    
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
    viewpointMid = floor((GenomicRanges::start(viewpointRegion) +
                          GenomicRanges::end(viewpointRegion)) / 2)
    
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
    for(i in seq_len(length(runTotalsList)))
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
    totals = GenomicRanges::mcols(runTotals)[,'runTotals']
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
    chr = as.character(GenomicRanges::seqnames(regionOfInterest)[1])
    fromHere = min(GenomicRanges::start(regionOfInterest))
    toHere = max(GenomicRanges::end(regionOfInterest))
    idx = (as.character(SummarizedExperiment::seqnames(deltaSE)) == chr &
           fromHere <= SummarizedExperiment::start(deltaSE) &
           SummarizedExperiment::end(deltaSE) <= toHere)
    deltaSE = deltaSE[idx,]

    ## ##########################################################################
    ## Bin to small bin size.  These give the start and end of each bin:
    start = seq(from=fromHere,to=toHere,by=smallBinSize)
    end = start - 1 + smallBinSize
    gr = GenomicRanges::GRanges(seqnames=chr,
                                IRanges::IRanges(start,end))

    
    binnedDeltaSE = binSummarizedExperiment(gr,deltaSE)

    ## ##########################################################################
    ## Bin to big bin size:
    lambda = bigBinSize / smallBinSize
    bigBinDeltaSE = rebinToMultiple(binnedDeltaSE,lambda)

    ## ##########################################################################
    ## Get real runs:
    realRunTotals = getRunTotals(bigBinDeltaSE,viewpointRegion)
    mid = floor((GenomicRanges::start(viewpointRegion) +
                 GenomicRanges::end(viewpointRegion)) / 2)
    realLopsidedness = getLopsidedness(bigBinDeltaSE,viewpointRegion,mid)

    scrambledDeltas = list()
    scrambledRuns = list()
    rowRanges=SummarizedExperiment::rowRanges(binnedDeltaSE)
    colData=SummarizedExperiment::colData(binnedDeltaSE)
    for(j in seq_len(numPermutations))
    {
        permutation = generatePermutation(binnedDeltaSE,viewpointRegion)
        m = matrix(SummarizedExperiment::assay(binnedDeltaSE)[permutation,],ncol=1)
        
        scrambledDeltaSE = SummarizedExperiment::SummarizedExperiment(colData=colData,
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
    answers = SummarizedExperiment::rowRanges(bigBinDeltaSE)

    mcols =
        data.frame(delta=SummarizedExperiment::assay(bigBinDeltaSE),
                   isViewpoint=FALSE,
                   lopsidednessIsSignificant=FALSE,
                   isSignificantByMin=FALSE,
                   isSignificantByMax=FALSE,
                   isSignificantByAbs=FALSE,
                   stringsAsFactors=FALSE)

    meets = IRanges::findOverlaps(answers,viewpointRegion)
    mcols$isViewpoint[IRanges::from(meets)] = TRUE

    ## ##########################################################################
    ## Is lopsidedness significant?
    if(realLopsidedness > cutoffs['lopsidedness'])
        mcols$lopsidednessIsSignificant[mcols$isViewpoint] = TRUE

    ## ##########################################################################
    ## What about runs?
    for(i in seq_len(length(realRunTotals)))
    {
        if(realRunTotals$runTotals[i] < cutoffs['min'])
        {
            meets = IRanges::findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByMin[IRanges::from(meets)] = TRUE
        }

        if(realRunTotals$runTotals[i] > cutoffs['max'])
        {
            meets = IRanges::findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByMax[IRanges::from(meets)] = TRUE
        }

        
        if(abs(realRunTotals$runTotals[i]) > cutoffs['abs'])
        {
            meets = IRanges::findOverlaps(answers,realRunTotals[i])
            mcols$isSignificantByAbs[IRanges::from(meets)] = TRUE
        }
    }
    GenomicRanges::mcols(answers) = mcols

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

    g = ggplot2::ggplot(df,ggplot2::aes(x=ctr,y=delta,fill=significance)) +
        ggplot2::geom_col() +
        ggplot2::ggtitle(title) +
        ggplot2::scale_discrete_manual(breaks=breaks,
                                       values=values,
                                       aesthetics=c('color','fill')) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=18),
                       axis.title=ggplot2::element_text(size=16,face="bold")) +
        ggplot2::xlab(xLabel)
    
    if(! legend)
        g = g + ggplot2::theme(legend.position='none')

    return(g)
}
