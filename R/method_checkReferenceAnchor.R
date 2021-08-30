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
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
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
