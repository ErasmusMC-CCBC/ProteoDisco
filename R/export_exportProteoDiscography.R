#' @title Export the generated (mutant) peptides sequences to a FASTA database.
#' @description Exports (mutant) peptide sequences to FASTA using several fields as identifier.
#' In addition, users can specify to only export (mutant) sequences containing at a min. number of proteotypic fragments.
#'
#' @param ProteoDiscography (\link[ProteoDisco]{ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param outFile (character): Filepath to output FASTA. If left NULL, it will return an AAStringSet with the given records.
#' @param minProteotypicFragments (integer): Only output mutant protein-isoforms from incorporated genomic variants with at least this many proteotypic fragments
#' (if `checkProteotypicFragments()` was performed).
#' @param aggregateSamples (logical): Should samples be aggregated into the same output database (TRUE) or should a seperate FASTA-file be generated per sample (FALSE).
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
#' ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
#'   ProteoDiscography = ProteoDiscography.hg19,
#'   # Do not aggregate samples and generate mutant transcripts from the mutations per sample.
#'   aggregateSamples = FALSE,
#'   # If there are multiple mutations within the same exon (CDS), place them on the same mutant CDS sequence.
#'   aggregateWithinExon = TRUE,
#'   # Aggregate multiple mutant exons (CDS) within the same transcripts instead of incorporating one at a time.
#'   aggregateWithinTranscript = TRUE,
#'   # If there are overlapping mutations on the same coding position, retain only the first of the overlapping mutations.
#'   # If set to FALSE, throw an error and specify which CDS had overlapping mutations.
#'   ignoreOverlappingMutations = TRUE,
#'   # Number of threads.
#'   threads = 1
#' )
#'
#' # Export peptide sequences to FASTA file. Optionally, only export those with at least a min.
#' # number of proteotypic fragments.
#' exportProteoDiscography(ProteoDiscography.hg19, outFile = 'out.fasta')
#'
#' @return Writes a FASTA file containing the mutant protein isoforms if outFile is given. Otherwise, will return an AAStringSet.
#' @importFrom rlang .data
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @export
exportProteoDiscography <- function(ProteoDiscography, outFile = NULL, minProteotypicFragments = NULL, aggregateSamples = TRUE) {

  # Input validation --------------------------------------------------------

  checkmate::assertClass(ProteoDiscography, classes = 'ProteoDiscography')
  checkmate::assertCharacter(outFile, null.ok = TRUE)
  checkmate::assertNumeric(minProteotypicFragments, null.ok = TRUE)
  checkmate::assertLogical(aggregateSamples)

  ParallelLogger::logTrace('ProteoDisco - Exporting mutant protein isoforms.')


  # Retrieve mutant transcript sequences and transform to tibble. -----------

  outputData <- ProteoDisco::mutantTranscripts(ProteoDiscography)
  outputData <- base::lapply(outputData, tibble::as_tibble)


  # Convert splice-junctions to one row per translated frame. ---------------

  if(base::length(outputData$spliceJunctions) != 0){

    outputData$spliceJunctions <- outputData$spliceJunctions %>%
      dplyr::select(c('identifierFASTA', 'AA.SequenceMut', 'sample')) %>%
      dplyr::group_split(.data$identifierFASTA)

    outputData$spliceJunctions <- dplyr::bind_rows(base::lapply(outputData$spliceJunctions, function(x){

      tibble::as_tibble(
        data.frame(
          identifierFASTA = x$identifierFASTA,
          sample = x$sample,
          frame = base::seq_len(S4Vectors::elementNROWS(x$AA.SequenceMut)),
          AA.SequenceMut = base::unlist(x$AA.SequenceMut)
        )
      ) %>% dplyr::mutate(identifier = sprintf('%s|Frame:%s', x$identifierFASTA, .data$frame))

    }))

  }

  # Only select sequences with (neo)proteotypic fragments -------------------

  # Check if checkProteotypicFragments was performed.
  if(!base::is.null(minProteotypicFragments)){

    checkGenomicVariants <- base::length(BiocGenerics::colnames(outputData$manualSequences)) != 0 & !'proteotypicFragmentsCount' %in% BiocGenerics::colnames(outputData$genomicVariants)
    checkManualSeqs <- base::length(BiocGenerics::colnames(outputData$manualSequences)) != 0 & !'proteotypicFragmentsCount' %in% BiocGenerics::colnames(outputData$manualSequences)

    if(base::any(checkGenomicVariants, checkManualSeqs)){

      stop('\tUser requested min. proteotypic fragments per output sequence but checkProteotypicFragments() was not yet performed.')

    }else{

      ParallelLogger::logInfo(sprintf('\tOnly exporting protein-isoforms with at least %s proteotypic fragment(s).', minProteotypicFragments))
      if(base::length(outputData$genomicVariants) != 0) {
        outputData$genomicVariants <- outputData$genomicVariants %>% dplyr::filter(.data$proteotypicFragmentsCount >= minProteotypicFragments)
      }
      if(base::length(outputData$manualSequences) != 0) {
        outputData$manualSequences <- outputData$manualSequences %>% dplyr::filter(.data$proteotypicFragmentsCount >= minProteotypicFragments)
      }
    }
  }


  # Aggregate output protein sequences --------------------------------------

  # Generate FASTA header (identifier) for the genomic mutations.
  outputData$genomicVariants <- outputData$genomicVariants %>%
    dplyr::mutate(
      identifier = base::sprintf(
        'MUT=%s;GENE=%s;TXID=%s;SAMPLE=%s',
        base::ifelse(.data$numberOfMutationsInTx > 2 | base::nchar(.data$mutations.genomic) >= 75,
                     paste(.data$numberOfMutationsInTx, 'mutation(s)'), .data$mutations.genomic),
        base::ifelse('SYMBOL' %in% colnames(.data), .data$SYMBOL, .data$geneID),
        .data$Tx.ID,
        .data$sample
      )
    )

  outputData$FASTA <- dplyr::bind_rows(
    outputData$genomicVariants %>% dplyr::select(c('identifier', 'AA.SequenceMut', 'sample'))
  )

  if(base::length(outputData$manualSequences) != 0){
    outputData$FASTA <- dplyr::bind_rows(
      outputData$FASTA,
      outputData$manualSequences %>% dplyr::select(c('identifier', 'AA.SequenceMut', 'sample'))
    )
  }

  if(base::length(outputData$spliceJunctions) != 0){
    outputData$FASTA <- dplyr::bind_rows(
      outputData$FASTA,
      outputData$spliceJunctions %>% dplyr::select(c('identifier', 'AA.SequenceMut', 'sample'))
    )
  }


  # Export to FASTA or AAStringSet ------------------------------------------

  ParallelLogger::logInfo(sprintf('\tExporting %s mutant protein sequences.', base::nrow(outputData$FASTA)))

  if(aggregateSamples){

    x <- Biostrings::AAStringSet(outputData$FASTA$AA.SequenceMut)
    base::names(x) <- outputData$FASTA$identifier

    if(!is.null(outFile)){

      ParallelLogger::logInfo(base::sprintf('\tWriting to %s.', outFile))
      Biostrings::writeXStringSet(x = x, filepath = outFile, append = FALSE, format = 'fasta')

    }

  } else {

    ParallelLogger::logInfo(sprintf('\t\tExporting into %s sample-specific fasta file(s).', dplyr::n_distinct(outputData$FASTA$sample)))

    # Loop per sample and write output into sample-specific FASTA with <sample>_<outFile> paths.
    x <- base::lapply(split(outputData$FASTA, outputData$FASTA$sample), function(s){

      x <- Biostrings::AAStringSet(s$AA.SequenceMut)
      names(x) <- s$identifier

      if(!is.null(outFile)){

        ParallelLogger::logInfo(base::sprintf('\tWriting to %s (per sample).', outFile))
        Biostrings::writeXStringSet(x = x, filepath = base::paste(unique(s$sample), outFile, sep = '_'), append = FALSE, format = 'fasta')

      }

      return(x)

    })

  }

  # Return statement --------------------------------------------------------

  base::ifelse(is.null(outFile), return(x), return(NULL))

}
