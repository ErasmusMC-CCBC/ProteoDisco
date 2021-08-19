require(GenomicFeatures)
require(VariantAnnotation)
require(BSgenome)
require(S4Vectors)
require(Biostrings)
require(tibble)

setClassUnion("DNAStringSet_OR_BSgenome", c("DNAStringSet", "BSgenome"))

# ProteoDiscography ----------------------------------------------------------------

#' ProteoDiscography
#'
#' @description An S4 object containing the reference genome sequences and gene-model annotations (TxDB).
#' 
#' This also stores genomic variants and splice-junctions from which mutant transcript sequences can be generated.
#'
#' @slot TxDb (\link[GenomicFeatures]{TxDb}): TxDb object containing the genomic and transcriptomic annotations.
#' @slot genomeSeqs (\link[Biostrings]{DNAStringSet}): Genomic sequence of the respective genome.
#' @slot input.genomicVariants (\link[VariantAnnotation]{VRangesList}): Imported genomic variants (SNV, MNV and InDels).
#' @slot input.spliceJunctions (DataFrame): Imported splice-junctions.
#' @slot input.manualSequences (DataFrame): Imported manual sequences.
#' @slot mutantTranscripts.genomicVariants (\link[tibble]{tibble}): Generated mutant mRNA sequences from genomic variants.
#' @slot mutantTranscripts.spliceJunctions (\link[tibble]{tibble}): Generated mutant mRNA sequences from splice-junctions.
#' @slot mutantTranscripts.manualSequences (\link[tibble]{tibble}): Processed mutant mRNA sequences from manual input.
#' @slot GENETIC_CODE (\link[Biostrings]{GENETIC_CODE_TABLE}): The genetic code table to be used during translation.
#' @slot metadata (data.frame): Supplied data.frame.
#' 
#' @rdname ProteoDiscography
#' @exportClass ProteoDiscography
#' @author Job van Riet
setClass(
  "ProteoDiscography", slots = methods::representation(
    
    # TxDb and genomic sequences.
    TxDb = 'TxDb',
    genomeSeqs = 'DNAStringSet_OR_BSgenome',
    
    # Stores the user-input.
    input.genomicVariants = 'VRangesList',
    input.spliceJunctions = 'tbl_df',
    input.manualSequences = 'DataFrame',
    
    # Stores the mutant transcripts we generated, checked and cleaned.
    # We stores these as DataFrames with sample and mutational information and a DNAString containing the mutant Tx sequence.
    mutantTranscripts.genomicVariants = 'DataFrame',
    mutantTranscripts.spliceJunctions = 'DataFrame',
    mutantTranscripts.manualSequences = 'DataFrame',
    GENETIC_CODE = 'character',
    
    # Stores additional information of the ProteoDiscography.
    metadata = 'data.frame'
  ),
  prototype = list(
    TxDb = NULL,
    genomeSeqs = NULL,
    input.genomicVariants = VariantAnnotation::VRangesList(),
    input.spliceJunctions  = tibble::tibble(),
    input.manualSequences  = S4Vectors::DataFrame(),
    mutantTranscripts.genomicVariants  = S4Vectors::DataFrame(),
    mutantTranscripts.spliceJunctions  = S4Vectors::DataFrame(),
    mutantTranscripts.manualSequences = S4Vectors::DataFrame(),
    GENETIC_CODE = Biostrings::GENETIC_CODE,
    metadata = data.frame('CreatedOn' = format(Sys.time(), "%a %b %d %H:%M:%S %Y"), ProteoDiscoVersion = packageVersion('ProteoDisco'), overlapUniqueCDS = NA, overlapUniqueGenes = NA)
  )
)
