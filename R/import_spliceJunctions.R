#' @title Import splice-junctions into the ProteoDiscograpy.
#' @description Generates putative gene-models based on supplied genomic coordinates of splice-junctions.
#'
#' Input should be a {tibble} containing the following columns:
##' \itemize{
##'  \item{junctionA: Genomic coordinates of the 5'-junction. (format: chr:start:strand, i.e.: chr1:100:+)}
##'  \item{junctionB: Genomic coordinates of the 3'-junction. (format: chr:end:strand, i.e.: chr1:150:+)}
##'  \item{sample: Names of the samples. (character, optional)}
##'  \item{identifier: The identifier which will be used in downstream analysis. (character, optional)}
##' }
#'
#' Common splice-junction formats (BED and SJ.out.tab (STAR)) can also be supplied and are converted into the correct DataFrame.
#'
#' By utilizing two separate junction-sites, interchromosomal trans-splicing or chimeric transcripts from genomic fusions
#' (e.g., resulting from the BCR/ABL1 fusion-gene) can also be handled.
#'
#' @param ProteoDiscography ({ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param inputSpliceJunctions (\link[tibble]{tibble}): Tibble containing the splice-junctions.
#' @param isTopHat (logical): Are the imported (.BED) files from TopHat? If so, the start-end of the SJ are corrected for max. overhang.
#' @param samples (character): Preferred names for the samples if BED / TAB files are supplied, default is derived from filepath.
#' @param aggregateSamples (logical): Should splice-junctions from multiple samples be aggregated? Or should sample-specific models be generated?
#' If genomic variants are to be incorporated within the derived splice-transcripts, the names of samples need to be match.
#' @param removeExisting (logical): Should existing entries be removed?
#' @param overwriteDuplicateSamples (logical): Should duplicate samples be overwritten?
#'
#' @examples
#'
#'  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#'    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'  )
#'
#'  # Import from file.
#'  ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
#'    ProteoDiscography = ProteoDiscography.hg19,
#'    inputSpliceJunctions = system.file('extdata', 'spliceJunctions_pyQUILTS_chr22.bed', package = 'ProteoDisco'),
#'    # (Optional) Rename samples.
#'    samples = 'pyQUILTS',
#'    # Specify that the given BED files are obtained from TopHat.
#'    # Chromosomal coordinates from TopHat require additional formatting.
#'    isTopHat = TRUE,
#'  )
#'
#'  # Or, import splice-junctions (even spanning different chromosomes) based on our format.
#'  testSJ <- readr::read_tsv(system.file('extdata', 'validationSetSJ_hg19.txt', package = 'ProteoDisco'))
#'
#'  # Add custom SJ to ProteoDiscography.
#'  ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
#'    ProteoDiscography = ProteoDiscography.hg19,
#'    inputSpliceJunctions = testSJ,
#'    # Append to existing SJ-input.
#'    removeExisting = FALSE
#'  )
#'
#' @return {ProteoDiscography} with imported splice-junctions.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#'
#' @export
importSpliceJunctions <- function(ProteoDiscography, inputSpliceJunctions, isTopHat = TRUE, samples = NULL, aggregateSamples = FALSE, removeExisting = FALSE, overwriteDuplicateSamples = FALSE){

  # Input validation. -------------------------------------------------------

  checkmate::assertClass(ProteoDiscography, classes = 'ProteoDiscography')
  checkmate::assertMultiClass(inputSpliceJunctions, classes = c('tbl_df', 'character'), null.ok = FALSE)
  checkmate::assertLogical(isTopHat)
  checkmate::assertCharacter(samples, null.ok = TRUE)
  checkmate::assertLogical(aggregateSamples)
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(overwriteDuplicateSamples)


  # Check the validity of the supplied tibble.
  if(!is.character(inputSpliceJunctions)){
    checkmate::assertCharacter(inputSpliceJunctions$junctionA)
    checkmate::assertCharacter(inputSpliceJunctions$junctionB)
    checkmate::assertCharacter(inputSpliceJunctions$sample, null.ok = TRUE)
    checkmate::assertCharacter(inputSpliceJunctions$identifier, null.ok = TRUE)
  }

  ParallelLogger::logInfo('ProteoDisco - Importing splice-junctions to the ProteoDiscography.')


  # Import splice-junctions -------------------------------------------------

  # Internal function for BED/TAB importation / conversion.
  importFile <- function(file, sample = NA, ProteoDiscography){

    data.File <- NULL

    # Import BED
    if(base::grepl('.bed$|.bed.gz$', file)){

      ParallelLogger::logTrace(sprintf('Importing BED (%s): %s', sample, file))

      data.File <- rtracklayer::import.bed(file)
      GenomeInfoDb::seqlevelsStyle(data.File) <- GenomeInfoDb::seqlevelsStyle(ProteoDiscography@TxDb)

      # Correct start / end for max. overhang in TopHat files.
      if(isTopHat){
        GenomicRanges::start(data.File) <- GenomicRanges::start(data.File) + base::unlist(base::lapply(GenomicRanges::width(data.File$blocks), function(x) x[1]))
        GenomicRanges::end(data.File) <- GenomicRanges::end(data.File) - base::unlist(base::lapply(GenomicRanges::width(data.File$blocks), function(x) x[2]))
      }

    }

    # Import SJ.out.tab files (STAR)
    if(base::grepl('.tab$', file)){

      ParallelLogger::logTrace(sprintf('Importing SJ.out.tab (%s): %s', sample, file))

      data.File <- readr::read_tsv(file, col_names = c('chr', 'start', 'end', 'strand', 'intronMotif', 'annotation', 'nUniqReads', 'nMultimapReads', 'maxOverhang'), col_types = 'ciiiiiiii')
      data.File <- data.File %>% dplyr::mutate(strand = dplyr::recode(.data$strand, '0' = '*', '1' = "+", '2' = "-"))
      data.File <- GenomicRanges::makeGRangesFromDataFrame(data.File, keep.extra.columns = TRUE)
      GenomeInfoDb::seqlevelsStyle(data.File) <- GenomeInfoDb::seqlevelsStyle(ProteoDiscography@TxDb)

    }

    if(!is.null(data.File)){

      # Convert to tibble.
      data.DF <- base::data.frame(
        sample = sample,
        identifier = NA,
        junctionA = base::sprintf('%s:%s:%s', as.character(GenomeInfoDb::seqnames(data.File)), GenomicRanges::start(data.File), GenomicRanges::strand(data.File)),
        junctionB = base::sprintf('%s:%s:%s', as.character(GenomeInfoDb::seqnames(data.File)), GenomicRanges::end(data.File), GenomicRanges::strand(data.File))
      ) %>% tibble::as_tibble()

      if(!is.null(data.File$name)) data.DF$identifier <- data.File$name

      return(data.DF)

    }else{

      ParallelLogger::logInfo(sprintf('The following file is not a .bed(.gz) or .tab file (%s): %s', sample, file))

      return(NULL)
    }

  }

  # Import files into DataFrame.
  if(is.character(inputSpliceJunctions)){
    data.Input <- dplyr::bind_rows(base::lapply(base::seq_len(base::length(inputSpliceJunctions)), function(i) importFile(file = inputSpliceJunctions[i], sample = samples[i], ProteoDiscography = ProteoDiscography)))
  }else{
    data.Input <- inputSpliceJunctions
  }

  # Aggregate samples.
  if(aggregateSamples){
    data.Input$sample <- 'Aggregated'
  }

  # Add to ProteoDiscography ------------------------------------------------

  # Append splice-junctions to the ProteoDiscography.
  ProteoDiscography <- .setSplicingJunctions(
    ProteoDiscography,
    data.Input,
    removeExisting,
    overwriteDuplicateSamples
  )

  # Return statement --------------------------------------------------------

  return(ProteoDiscography)
}
