#' @title Checks if mutations and reference genome correspond
#' @description Overlaps the reference bases and positions with the given genome and check if all positions overlap and reports
#' non-overlapping mutations.
#' 
#' @param mutations (\link[VariantAnnotation]{VRanges}): Mutations to incorporate within transcripts.
#' @param genomeSeqs (\link[BSgenome]{BSgenome} or \link[Biostrings]{DNAStringSet}): Reference genome sequences.
#' @param ignoreNonMatch (logical): Should non-matching reference anchors be ignored? These mutations will be removed.
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
#' @return (logical): Returns TRUE/FALSE if per given mutation based on matching position and base with the reference genome.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
.checkReferenceAnchor <- function(mutations, genomeSeqs, ignoreNonMatch = FALSE, threads = 1) {
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertClass(mutations, classes = 'VRanges')
  checkmate::assert(
    checkmate::checkClass(genomeSeqs, classes = 'DNAStringSet'),
    checkmate::checkClass(genomeSeqs, classes = 'BSgenome'), combine = 'or'
  )
  checkmate::checkInt(threads)
  
  ParallelLogger::logTrace(sprintf('ProteoDisco - Checking sanity of reference anchors for %s.', unique(mutations$sample)))
  
  # Overlap the reference anchors -------------------------------------------
  
  # Per chromosome, check the reference position and base.
  checkRef <- base::do.call(base::rbind, base::lapply(GenomeInfoDb::seqlevelsInUse(mutations), function(chr){
    
    ParallelLogger::logTrace(sprintf('ProteoDisco - Checking sanity of reference anchors in %s (using %s threads)', chr, threads))
    
    g <- genomeSeqs[[chr]]
    mutsChr <- mutations[GenomeInfoDb::seqnames(mutations) == chr,]
    
    base::do.call(base::rbind, BiocParallel::bplapply(base::seq_along(mutsChr), function(i){
      mut <- mutsChr[i]
      base::data.frame(
        refBase = base::as.character(g[IRanges::start(mut):IRanges::end(mut)]),
        refCheck = g[IRanges::start(mut):IRanges::end(mut)] == Biostrings::DNAString(VariantAnnotation::ref(mut))
      )
    }, BPPARAM = BiocParallel::MulticoreParam(workers = threads)))
  }))
  
  # Check if any reference anchor did not correspond to the reference genome.
  if(!all(checkRef$refCheck)){
    notCorrect <- mutations[!checkRef$refCheck]
    
    message <- base::sprintf(
      '%s supplied reference bases and positions did not match with the given reference genome (%s) for %s. %s', 
      base::length(notCorrect), 
      base::unique(GenomeInfoDb::genome(genomeSeqs)), 
      base::unique(notCorrect$sample),
      base::ifelse(ignoreNonMatch, 'Will remove these records and continue.', 'You can remove these elements with ignoreNonMatch = TRUE to continue but it is advised to check these positions.'))
    
    if(ignoreNonMatch) warning(message)
    if(!ignoreNonMatch) stop(message)
  }
  
  # Return a list of TRUE/FALSE based on matching reference anchor with given genome.
  return(checkRef$refCheck)
  
}
