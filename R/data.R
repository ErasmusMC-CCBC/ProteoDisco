#' Example ProteoDiscography.
#'
#' A ProteoDiscography object based on validated genomic variants (hg19) using existing transcript annotations (hg19).
#'
#' Generated using the following commands:
#' # Use existing hg19 annotations.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#'
#' # Add mutations from hg19 validation set.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::importGenomicVariants(
#'   ProteoDiscography = ProteoDiscographyExample.hg19,
#'   files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'),
#'   samplenames = 'Validation Set (GRCh37)',
#'   threads = 1
#' )
#'
#' # Incorporate genomic mutations.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::incorporateGenomicVariants(
#'   ProteoDiscography = ProteoDiscographyExample.hg19,
#'   aggregateSamples = FALSE,
#'   aggregateWithinExon = TRUE,
#'   aggregateWithinTranscript = FALSE,
#'   ignoreOverlappingMutations = TRUE,
#'   threads = 1
#' )
#'
#' @source \url{https://cancer.sanger.ac.uk/cosmic}
"ProteoDiscographyExample.hg19"
