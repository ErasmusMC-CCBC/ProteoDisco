#' @title Generate the annotation database (ProteoDiscography)
#' @description Generates a database containing the genomic and transcriptomic annotations which will be used downstream 
#' to integrate genomic variants and discover mutant protein sequences with ProteoDisco.
#' 
#' By default, only the common seqlevels between the TxDb and genomeSeqs will be kept..
#'
#' @param TxDb (\link[GenomicFeatures]{TxDb}): TxDb object containing the genomic and transcriptomic annotations.
#' @param genomeSeqs (\link[Biostrings]{DNAStringSet} or \link[BSgenome]{BSgenome}): Genomic sequence of the respective genome.
#' @param useOnlySharedSeqlevels (logical): Should only shared seqlevels between the TxDb and genomeSeqs be kept?
#' @param geneticCode (character): Which codon-table should be used? Default is set to Standard, check \link[Biostrings]{GENETIC_CODE_TABLE} for additional options.
#' 
#' @examples
#'  # Generate a ProteoDiscography using existing TxDb and annotations.
#'  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'  )
#'  
#' @return ProteoDiscography
#' @concept methods
#' @keywords methods
#' @importFrom methods new
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @export
generateProteoDiscography <- function(TxDb, genomeSeqs, useOnlySharedSeqlevels = TRUE, geneticCode = 'Standard'){
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertClass(TxDb, classes = 'TxDb')
  checkmate::assert(
    checkmate::checkClass(genomeSeqs, classes = 'DNAStringSet'),
    checkmate::checkClass(genomeSeqs, classes = 'BSgenome'), combine = 'or'
  )
  checkmate::assertCharacter(geneticCode)
  
  # Check if genetic code exists.
  checkmate::assertCharacter(Biostrings::getGeneticCode(id_or_name2 = geneticCode, full.search = TRUE))
  
  # Retrieve common seqlevels and subset data accordingly.
  if(useOnlySharedSeqlevels & methods::is(genomeSeqs, 'DNAStringSet')){
    ParallelLogger::logTrace('ProteoDisco - Subsetting on shared seqlevels.')
    
    genomeSeqs <- genomeSeqs[names(genomeSeqs) %in% GenomeInfoDb::seqlevels(TxDb)]
    TxDb <- GenomeInfoDb::keepSeqlevels(TxDb, names(genomeSeqs), pruning.mode = 'coarse')
    
  }
  
  if(useOnlySharedSeqlevels & methods::is(genomeSeqs, 'BSgenome')) ParallelLogger::logTrace('ProteoDisco - Ignoring subsetting on shared seqlevels as a BSgenome was given.')
  
  ParallelLogger::logInfo('ProteoDisco - Generating the annotation database (ProteoDiscography)')
  
  
  # Sanity-checks -----------------------------------------------------------
  
  # Check for common seqlevels; ignore if a BSgenome was given.
  if(methods::is(genomeSeqs, 'DNAStringSet') & useOnlySharedSeqlevels){
    if(useOnlySharedSeqlevels & !all(GenomeInfoDb::seqlevels(TxDb) %in% names(genomeSeqs))) stop(sprintf('The given TxDb and genomeSeqs have incompatible seqlevels: %s', paste(GenomeInfoDb::seqlevels(TxDb)[!GenomeInfoDb::seqlevels(TxDb) %in% names(genomeSeqs)], collapse = ', ')))
    if(!useOnlySharedSeqlevels & !all(GenomeInfoDb::seqlevels(TxDb) %in% names(genomeSeqs))) warning(sprintf('The given TxDb and genomeSeqs have incompatible seqlevels but this will be ignored: %s', paste(GenomeInfoDb::seqlevels(TxDb)[!GenomeInfoDb::seqlevels(TxDb) %in% names(genomeSeqs)], collapse = ', ')))
  }
  
  # Create new ProteoDiscography -----------------------------------------------------
  
  x <- new('ProteoDiscography', TxDb = TxDb, genomeSeqs = genomeSeqs)
  
  # If requested, change default genetic code table.
  if(geneticCode != 'Standard') x@GENETIC_CODE <- Biostrings::getGeneticCode(id_or_name2 = geneticCode, full.search = TRUE)
  
  return(x)
}

