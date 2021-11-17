# Internal setter methods.
setGeneric('.setGenomicVariants', function(x, variants, removeExisting, overwriteDuplicateSamples) standardGeneric(".setGenomicVariants"))
setGeneric('.setManualSequences', function(x, transcripts, removeExisting, overwriteDuplicateSamples) standardGeneric(".setManualSequences"))
setGeneric('.setSplicingJunctions', function(x, spliceJunctions, removeExisting, overwriteDuplicateSamples) standardGeneric(".setSplicingJunctions"))

#' @title Change the underlying TxDb of a ProteoDiscography object.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @param TxDb (\link[GenomicFeatures]{TxDb}): TxDb object containing the genomic and transcript annotations.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' summary(ProteoDiscographyExample.hg19)
#'
#' @return ProteoDiscography with updated TxDb.
#' @rdname setTxDb
#' @export
setGeneric("setTxDb", function(x, TxDb) standardGeneric("setTxDb"))

#' @title Change the underlying genomic sequences of a ProteoDiscography object.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @param genomeSeqs (\link[Biostrings]{DNAStringSet} or \link[BSgenome]{BSgenome}): Genomic sequence of the respective genome.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link genomic sequences.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' summary(ProteoDiscographyExample.hg19)
#'
#' @return ProteoDiscography with updated genomic sequences.
#' @rdname setGenomicSequences
#' @export
setGeneric("setGenomicSequences", function(x, genomeSeqs) standardGeneric("setGenomicSequences"))

#' @title Retrieve imported genomic variants, splice-junctions and manual sequences.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Retrieve the imported records.
#' getDiscography(ProteoDiscographyExample.hg19)
#'
#' @return Return a list of imported records, per category (genomic variants, splice-junctions and manual sequences).
#' @rdname getDiscography
#' @export
setGeneric('getDiscography', function(x) standardGeneric("getDiscography"))

#' @title Adds mutant transcript sequences to the ProteoDiscography in the appropriate slot
#' @param x (ProteoDiscography): The ProteoDiscography for which the slot will be edited.
#' @param transcripts (DataFrame): Transcripts to be used in the slot.
#' @param slotType (character): Implemented slot to be edited.
#'
#' @examples
#'
#' # From a ProteoDiscography with imported and incorporated records, take only the first 10 records.
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography)$genomicVariants[1:10], slotType = 'genomicVariants')
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography)$spliceJunctions[1:10], slotType = 'spliceJunctions')
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography)$manualSequences[1:10], slotType = 'manualSequences')
#'
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Only keep the first ten records.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::setMutantTranscripts(
#'   x = ProteoDiscographyExample.hg19,
#'   transcripts = ProteoDisco::mutantTranscripts(ProteoDiscographyExample.hg19)$genomicVariants[1:10,],
#'   slotType = 'genomicVariants'
#' )
#'
#' @return {ProteoDiscography} with updated records.
#' @rdname setMutantTranscripts
#' @export
setGeneric("setMutantTranscripts", function(x, transcripts, slotType) standardGeneric("setMutantTranscripts"))

#' @title Adds mutant transcript sequences to the ProteoDiscography in the appropriate slot
#' @param x (ProteoDiscography): ProteoDiscography.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Retrieve all generated mutant transcripts.
#' mutantTranscripts(ProteoDiscographyExample.hg19)
#'
#' @return Return all incorporated mutant transcripts (list of DataFrame).
#' @rdname mutantTranscripts
#' @export
setGeneric("mutantTranscripts", function(x) standardGeneric("mutantTranscripts"))

#' @title Determine proteotypic fragments within the translated CDS
#' @description Proteotypic fragments are checked against the input TxDb (and additional peptide sequences) to infer whether peptide sequences
#' contain one or multiple proteotypic fragment(s) after cleavage by a selected protease.
#'
#' @param x ({ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param enzymUsed (character): Preferred proteasome used in cleaving the peptide sequences, e.g. trypsin.
#' @param missedCleavages (integer): Number of subsequent missed cleavages.
#' @param additionalPeptides (AAStringSet): Additional peptide sequences against which to check proteotypic fragments, besides the input TxDb.
#' If left empty, only the (translated) input TxDb will be checked.
#' @param checkWithinMutantSeqs (logical): Should proteotypic fragments also not be present within another mutant peptide-sequence?
#'
#' @return ProteoDiscography with additional information added to the transcript sequence information.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Results will now contain additional information about proteotypic fragments.
#' # With a large TxDb, this can take a while.
#' # ProteoDiscography.hg19 <- ProteoDisco::checkProteotypicFragments(ProteoDiscographyExample.hg19)
#'
#' @return {ProteoDiscography} with an additional column specifying the number of proteotypic fragments per record.
#' @rdname checkProteotypicFragments
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @import cleaver
#' @export
setGeneric("checkProteotypicFragments", function(x, enzymUsed = 'trypsin', missedCleavages = 0, additionalPeptides = NULL, checkWithinMutantSeqs = FALSE) standardGeneric("checkProteotypicFragments"))
