#' @title Incorporate manual transcript sequences into the ProteoDiscography
#' @description The given DataFrame (transcripts) should include the following columns:
##' \itemize{
##'  \item{sample: Names of the samples. (character)}
##'  \item{identifier: The identifier which will be used in downstream analysis. (character)}
##'  \item{Tx.SequenceMut: The mutant mRNA sequence. (DNAStringSet)}
##'  \item{gene: Name of the gene. (character, optional)}
##' } 
##' 
#' The supplied mutant sequences do not need to have a respective record within the used TxDb as these are handled as independent sequences.
#' 
#' @inheritParams importGenomicVariants
#' @param transcripts (\link[S4Vectors]{DataFrame}): DataFrame containing metadata and the transcript sequences.
#' 
#' @examples
#' 
#'  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'  )
#'  
#'  manualSeq <- S4Vectors::DataFrame(
#'    sample = rep('Example', 2),
#'    identifier = c('Example 1', 'Example 2'),
#'    gene = c('GeneA', 'GeneB'),
#'    Tx.SequenceMut = Biostrings::DNAStringSet(list(Biostrings::DNAString('ATCGGGCCCGACGTT'), Biostrings::DNAString('GCTAGCGATCAGGGA')))
#'  )
#' 
#'  # Add to ProteoDiscography.
#'  ProteoDiscography.hg19 <- ProteoDisco::importTranscriptSequences(
#'    ProteoDiscography.hg19, 
#'    transcripts = manualSeq
#'  )
#'
#' @return {ProteoDiscography} with added mutant transcript sequences from manual input.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' 
#' @export
importTranscriptSequences <- function(ProteoDiscography, transcripts, removeExisting = FALSE, overwriteDuplicateSamples = FALSE){
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertClass(ProteoDiscography, classes = "ProteoDiscography")
  checkmate::assertClass(transcripts, classes = 'DataFrame')
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(overwriteDuplicateSamples)
  
  # Check if the given columns are of the correct format.
  checkmate::assertCharacter(transcripts$sample)
  checkmate::assertCharacter(transcripts$identifier)
  checkmate::assertCharacter(transcripts$gene, null.ok = TRUE)
  checkmate::assertClass(transcripts$Tx.SequenceMut, classes = 'DNAStringSet')
  
  ParallelLogger::logInfo('ProteoDisco - Importing manual transcript sequences to the ProteoDiscography.')
  
  
  # Import transcript sequences ---------------------------------------------
  
  # Append sequences to the ProteoDiscography.
  ProteoDiscography <- .setManualSequences(
    ProteoDiscography, 
    transcripts,
    removeExisting, 
    overwriteDuplicateSamples
  )
  
  
  # Return statement --------------------------------------------------------
  
  return(ProteoDiscography)
}
