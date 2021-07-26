#' @title Import genomic variants into the ProteoDiscography
#' @description Imports genomic variants (SNV, MNV and InDels) present within the supplied VCF/MAF files 
#' into the {ProteoDiscography} as a {VRanges}. This genomic variants can later be incorporated within transcript 
#' sequences at a later stage.
#'
#' @param ProteoDiscography ({ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param files (character): Path(s) to VCF or MAF files.
#' @param samplenames (character): Descriptive samplename(s) of the VCF files in the same order as input VCF file(s), if NULL the basename of the file will be used instead.
#' @param removeExisting (logical): Should previous mutations within the ProteoDiscography be removed?
#' @param overwriteDuplicateSamples (logical): Replace duplicate samples (TRUE) or throw an error if duplicate samples are found.
#' @param performAnchorCheck (logical): Should the reference anchor be check for consistency with the given genomic sequences?
#' @param ignoreNonMatch (logical): Should non-matching reference anchors be ignored? These mutations will be removed prior to appending.
#' @param threads (integer): Number of threads.
#' 
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' # Supply the ProteoDiscography with genomic variants to incorporate in downstream analysis. This can be one or multiple VCF / MAF files.
#' # Additional manual sequences and exon-exon mapping (i.e., splice junctions) can also be given as shown in the sections below.
#' ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
#'   ProteoDiscography = ProteoDiscography.hg19,
#'   # Provide the VCF / MAF files, if more then one supply a vector of files and corresponding samplenames.
#'   files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'), 
#'   # We can replace the original samples within the VCF with nicer names.
#'   samplenames = 'Validation Set (GRCh37)',
#'   # Number of threads used for parallelization.
#'   # We run samples sequentially and parallelize within (variant-wise multi-threading).
#'   threads = 1, 
#'   # To increase import-speed for this example, do not check for validity of the reference anchor with the given reference sequences.
#'   performAnchorCheck = FALSE
#' )
#'  
#' @return {ProteoDiscography} with additional imported SNVs, MNVs and InDels.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @export
importGenomicVariants <- function(ProteoDiscography, files, samplenames = NULL, removeExisting = FALSE, overwriteDuplicateSamples = TRUE, performAnchorCheck = TRUE, ignoreNonMatch = FALSE, threads = 1) {
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertFile(files, access = 'r')
  checkmate::assertCharacter(samplenames, len = length(files), null.ok = TRUE)
  checkmate::assertLogical(overwriteDuplicateSamples)
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(performAnchorCheck)
  checkmate::assertLogical(ignoreNonMatch)
  checkmate::assertInt(threads)
  
  if( ! base::all(base::grepl('\\.vcf$|\\.vcf.gz$|\\.maf$', files)) ) stop('Not all given files are either VCF, VCF.gz or MAF files!')
  
  # If no samplenames were specified, use the filenames as samplenames. 
  if(is.null(samplenames)) samplenames <- base::basename(files)
  
  ParallelLogger::logInfo(sprintf('ProteoDisco - Importing genomic variants from %d file(s) to the ProteoDiscography.', base::length(files)))
  
  
  # Import the mutations based on the file extension ------------------------
  
  # Initialize placeholders.  
  muts.VCF <- VariantAnnotation::VRangesList(); muts.MAF <- VariantAnnotation::VRangesList()
  
  if(any(grepl('\\.vcf$|\\.vcf.gz$', files))) muts.VCF <- importGenomicVariants.VCF(ProteoDiscography, files[grepl('\\.vcf$|\\.vcf.gz$', files)], samplenames[grepl('\\.vcf$|\\.vcf.gz$', files)], performAnchorCheck = performAnchorCheck, ignoreNonMatch = ignoreNonMatch, threads = threads)
  if(any(grepl('\\.maf$|\\.maf.gz$', files))) muts.MAF <- importGenomicVariants.MAF(ProteoDiscography, files[grepl('\\.maf$|\\.maf.gz$', files)], samplenames[grepl('\\.maf$|\\.maf.gz$', files)], performAnchorCheck = performAnchorCheck, ignoreNonMatch = ignoreNonMatch, threads = threads)
  
  # Append new samples (with genomic mutations) to the ProteoDiscography, append new mutations to existing samples if overwriteDuplicateSamples = FALSE.
  ProteoDiscography <- .setGenomicVariants(
    ProteoDiscography, 
    variants = c(muts.VCF, muts.MAF), 
    removeExisting = removeExisting, 
    overwriteDuplicateSamples = overwriteDuplicateSamples
  )
  
  # Calculate (or update) a quick overlap of genomic mutations with unique CDS and genes.
  ParallelLogger::logInfo('\tProteoDisco - Calculating a prelim. overlap of genomic mutations with given annotations')
  
  allMuts <- base::unique(GenomicRanges::GRanges(base::unlist(ProteoDiscography@input.genomicVariants)))
  
  allMutsOverlap <- GenomicFeatures::cdsByOverlaps(ProteoDiscography@TxDb, allMuts, minoverlap = 1, type = 'any', column = c('CDSID', 'GENEID'))
  ProteoDiscography@metadata$overlapUniqueCDS <- base::length(base::unique(allMutsOverlap$CDSID))
  ProteoDiscography@metadata$overlapUniqueGenes <- base::length(base::unique(unlist(allMutsOverlap$GENEID)))
  
  
  # Return statement --------------------------------------------------------
  
  return(ProteoDiscography)
  
}


# Internal function to check validity of import ---------------------------

.checkValidity.Import <- function(ProteoDiscography, data, samplename, performAnchorCheck, ignoreNonMatch, threads){
  
  # Use standarized chromosomal names.
  GenomeInfoDb::seqlevelsStyle(data) <- 'UCSC'
  
  # Fix weird genome names.
  GenomeInfoDb::genome(data) <- gsub('\"', '', GenomeInfoDb::genome(data))
  
  # Check if all mutations lie upon given chromosomes.
  checkChromosomes <- GenomeInfoDb::seqlevelsInUse(data) %in% GenomeInfoDb::seqlevels(ProteoDiscography@genomeSeqs)
  if(!all(checkChromosomes)) stop('There are mutations on chromosomes which were not given as genomic sequences: ', base::paste(base::unique(base::as.character(GenomeInfoDb::seqnames(data)[!checkChromosomes])), collapse = ', '))
  
  # Check correct REF basesas we use these as later anchors.
  if(any(!grepl('[ATCGN]', VariantAnnotation::ref(data)))) stop(base::sprintf('There are non-[ATCGN] bases in your reference nucleotides within file: %s', file))
  
  # Determine mutational type.
  data$mutationalType <- 'Other'
  if(any(VariantAnnotation::isSNV(data))) data[VariantAnnotation::isSNV(data),]$mutationalType <- 'SNV'
  if(any(VariantAnnotation::isIndel(data))) data[VariantAnnotation::isIndel(data),]$mutationalType <- 'InDel'
  if(any(!VariantAnnotation::isSNV(data) & !VariantAnnotation::isIndel(data))) data[(!VariantAnnotation::isSNV(data) & !VariantAnnotation::isIndel(data)),]$mutationalType <- 'MNV'
  data$mutationalType <- factor(data$mutationalType)
  
  # Add the sample name to the data.
  data$sample <- base::factor(samplename)
  
  # Fix sorting.
  data <- GenomeInfoDb::sortSeqlevels(data)
  data <- BiocGenerics::sort(data)
  
  if(performAnchorCheck){
    # Check reference anchors and remove non-matching reference anchors.
    correctAnchor <- .checkReferenceAnchor(data, ProteoDiscography@genomeSeqs, ignoreNonMatch = ignoreNonMatch, threads = threads)
    data <- data[correctAnchor]
  }else{
    ParallelLogger::logTrace(sprintf('ProteoDisco - Skipping reference-anchor checking for: %s', samplename))
  }
  
  # Return statement.
  return(data)
  
}

#' @title Import VCF files into a ProteoDiscography.
#' @inheritParams importGenomicVariants
#' 
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' # Supply the ProteoDiscography with genomic variants to incorporate in downstream analysis. This can be one or multiple VCF / MAF files.
#' # Additional manual sequences and exon-exon mapping (i.e., splice junctions) can also be given as shown in the sections below.
#' ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
#'   ProteoDiscography = ProteoDiscography.hg19,
#'   # Provide the VCF / MAF files, if more then one supply a vector of files and corresponding samplenames.
#'   files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'), 
#'   # We can replace the original samples within the VCF with nicer names.
#'   samplenames = 'Validation Set (GRCh37)',
#'   # Number of threads used for parallelization.
#'   # We run samples sequentially and parallelize within (variant-wise multi-threading).
#'   threads = 1, 
#'   # To increase import-speed for this example, do not check for validity of the reference anchor with the given reference sequences.
#'   performAnchorCheck = FALSE
#' )
#'  
#' @return {ProteoDiscography} with additional imported SNVs, MNVs and InDels.
#' @export
importGenomicVariants.VCF <- function(ProteoDiscography, files, samplenames = NULL, performAnchorCheck = TRUE, ignoreNonMatch = FALSE, threads = 1) {
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertFile(files, access = 'r')
  checkmate::assertCharacter(samplenames, len = length(files), null.ok = TRUE)
  checkmate::assertLogical(performAnchorCheck)
  checkmate::assertLogical(ignoreNonMatch)
  checkmate::assertInt(threads)
  
  if( ! all(grepl('\\.vcf|\\.vcf.gz', files)) ) stop('Not all given files are VCF or VCF.gz files!')
  
  # If no samplenames were specified, use the filenames as samplenames. 
  if(is.null(samplenames)) samplenames <- base::make.unique(base::basename(files))
  
  ParallelLogger::logInfo(sprintf('ProteoDisco - Importing %d VCF file(s).', base::length(files)))
  
  
  # Import VCF data to GRangesList ------------------------------------------
  
  # Read VCF file without INFO, GENO and SAMPLE fields to greatly speed up reading.
  param <- VariantAnnotation::ScanVcfParam(info = NA, geno = NA)
  
  # Read the VCF files into a VRangesList.
  data.VCF <- base::suppressWarnings(VariantAnnotation::VRangesList(
    BiocParallel::bplapply(files, function(file){
      
      ParallelLogger::logTrace(sprintf('ProteoDisco - Importing %s', file))
      
      # Import
      data <- VariantAnnotation::readVcfAsVRanges(file, param = param)
      
      # Check validity.
      data <- .checkValidity.Import(ProteoDiscography, data, samplenames[base::min(base::which(files == file))], performAnchorCheck, ignoreNonMatch, threads)
      
      # Return validated data.
      return(data)
      
    }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, stop.on.error = TRUE, progressbar = TRUE))
  ))
  
  # Add samplenames which will later be used to track from which file/sample a mutation originated.
  if(!is.null(samplenames)){
    base::names(data.VCF) <- samplenames
  }else{
    base::names(data.VCF) <- basename(files)
  }
  
  
  # Return statement --------------------------------------------------------
  
  return(data.VCF)
  
}


#' @title Import MAF files into a ProteoDiscography.
#' @inheritParams importGenomicVariants
#' 
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' # Supply the ProteoDiscography with genomic variants to incorporate in downstream analysis. This can be one or multiple VCF / MAF files.
#' # Additional manual sequences and exon-exon mapping (i.e., splice junctions) can also be given as shown in the sections below.
#' ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
#'   ProteoDiscography = ProteoDiscography.hg19,
#'   # Provide the VCF / MAF files, if more then one supply a vector of files and corresponding samplenames.
#'   files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'), 
#'   # We can replace the original samples within the VCF with nicer names.
#'   samplenames = 'Validation Set (GRCh37)',
#'   # Number of threads used for parallelization.
#'   # We run samples sequentially and parallelize within (variant-wise multi-threading).
#'   threads = 1, 
#'   # To increase import-speed for this example, do not check for validity of the reference anchor with the given reference sequences.
#'   performAnchorCheck = FALSE
#' )
#'  
#' @return {ProteoDiscography} with additional imported SNVs, MNVs and InDels.
#' @export
importGenomicVariants.MAF <- function(ProteoDiscography, files,  performAnchorCheck = TRUE, ignoreNonMatch = FALSE, samplenames = NULL, threads = 1) {
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertFile(files, access = 'r')
  checkmate::assertCharacter(samplenames, len = length(files), null.ok = TRUE)
  checkmate::assertLogical(performAnchorCheck)
  checkmate::assertLogical(ignoreNonMatch)
  checkmate::assertInt(threads)
  
  if( ! all(grepl('\\.maf|\\.maf.gz', files)) ) stop('Not all given files are MAF or MAF.gz files!')
  
  # If no samplenames were specified, use the filenames as samplenames. 
  if(is.null(samplenames)) samplenames <- base::make.unique(base::basename(files))
  
  ParallelLogger::logInfo(sprintf('ProteoDisco - Importing %d MAF file(s).',  base::length(files)))
  
  
  # Import MAF data to GRangesList ------------------------------------------ 
  
  data.MAF <- BiocParallel::bplapply(files, function(file){
    
    # Import raw data.
    data <- readr::read_delim(file, delim = '\t')
    
    # Convert to GRanges.
    data.GRanges <- GenomicRanges::makeGRangesFromDataFrame(data, start.field = 'Start_Position', end.field = 'End_Position', keep.extra.columns = TRUE)
    
    # Convert to VRanges.
    data.VRanges <- VariantAnnotation::makeVRangesFromGRanges(data.GRanges, ref.field = 'Reference_Allele', alt.field = 'Tumor_Seq_Allele2')
    
    # Check validity.
    data <- .checkValidity.Import(ProteoDiscography, data.VRanges, samplenames[base::min(base::which(files == file))], performAnchorCheck, ignoreNonMatch, threads)
    
    # Return validated data.
    return(data)
    
  }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, stop.on.error = TRUE, progressbar = TRUE))
  
  # Add samplenames which will later be used to track from which file/sample a mutation originated.
  if(!is.null(samplenames)){
    base::names(data.MAF) <- samplenames
  }else{
    base::names(data.MAF) <- base::basename(files)
  }
  
  
  # Return statement --------------------------------------------------------
  
  return(data.MAF)
}
