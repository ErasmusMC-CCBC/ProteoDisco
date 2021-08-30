#' @title Display a summary of ProteoDiscography
#' @param object (ProteoDiscography): ProteoDiscography object.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' print(ProteoDiscographyExample.hg19)
#'
#' @return Displays information relating to the ProteoDiscography into the console.
#' @rdname print
#' @export
setGeneric('print', function(object) standardGeneric("print"))

# Internal setter methods.
setGeneric('.setGenomicVariants', function(x, variants, removeExisting, overwriteDuplicateSamples) standardGeneric(".setGenomicVariants"))
setGeneric('.setManualSequences', function(x, transcripts, removeExisting, overwriteDuplicateSamples) standardGeneric(".setManualSequences"))
setGeneric('.setSplicingJunctions', function(x, spliceJunctions, removeExisting, overwriteDuplicateSamples) standardGeneric(".setSplicingJunctions"))

#' @title Display a summary of ProteoDiscography
#' @param object (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' show(ProteoDiscographyExample.hg19)
#'
#' @return Displays information relating to the ProteoDiscography into the console.
#' @rdname show
#' @export
setGeneric('show', function(object) standardGeneric("show"))

#' @title Display a summary of ProteoDiscography
#' @param object (ProteoDiscography): ProteoDiscography object.
#' @param verbose (logical): Set the verbosity
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' summary(ProteoDiscographyExample.hg19)
#'
#' @return Displays information relating to the ProteoDiscography into the console.
#' If verbose is set to TRUE, it will also output the imported genomic variants, splicing junctions and manual sequences as a list of tibbles.
#' @rdname summary
#' @export summary
setGeneric('summary', function(object, verbose = TRUE) standardGeneric("summary"))

#' @title Retrieve seqinfo.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' seqinfo(ProteoDiscographyExample.hg19)
#'
#' @return Returns an Seqinfo-class containing information on the imported genomic sequences.
#' @rdname seqinfo
#' @export
setGeneric('seqinfo', function(x) standardGeneric("seqinfo"))

#' @title Retrieve seqlevels.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' seqlevels(ProteoDiscographyExample.hg19)
#'
#' @return seqlevels of the improted genomic sequences.
#' @rdname seqlevels
#' @export
setGeneric('seqlevels', function(x) standardGeneric("seqlevels"))

#' @title Display the organism of the TxDb
#' @param object (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' organism(ProteoDiscographyExample.hg19)
#'
#' @return Returns the specified organism during creation.
#' @rdname organism
#' @export
setGeneric('organism', function(object) standardGeneric("organism"))

#' @title Retrieve imported genomic variants, splice-junctions and manual sequences.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @examples
#'
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
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
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' # Only keep the first ten records.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::setMutantTranscripts(
#'   x = ProteoDiscography.hg19,
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
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
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
#' # Import example ProteoDiscography (hg19) and re-link TxDb.
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
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
