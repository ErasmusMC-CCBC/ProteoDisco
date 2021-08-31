#' @title Incorporate genomic events into their overlapping exonic sequences
#'
#' @description \code{incorporateMutations} Incorporates SNV, MNV and InDels present in the {ProteoDiscography} on the transcripts.
#'
#' @param ProteoDiscography ({ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param aggregateSamples (logical): Should genomic mutations from different samples be incorporated within the same transcript?
#' @param aggregateWithinExon (logical): Should multiple mutations within the same exon be aggregated (TRUE) or should each mutation per exon produce a separate mutant transcript?
#' @param aggregateWithinTranscript (logical): Should multiple mutant exons within the same transcript be aggregated?
#' @param ignoreOverlappingMutations (logical): Incorporate first mutation (by order) and ignore subsequent overlapping mutations (and provide a warning) or stop the incorporation (if set to TRUE).
#' If aggregateWithinExon is set to FALSE, each mutations are added into separate exon to produce seperate transcripts.
#' @param threads (integer): Number of threads.
#'
#' @examples
#'
#' # Import example ProteoDiscography (hg19)
#' data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
#' ProteoDiscographyExample.hg19 <- setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ProteoDiscographyExample.hg19 <- setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Incorporate the genomic variants.
#' ProteoDiscographyExample.hg19 <- ProteoDisco::incorporateGenomicVariants(
#'   ProteoDiscography = ProteoDiscographyExample.hg19,
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
#' @return {ProteoDiscography} with mutant transcript sequences containing SNVs, MNVs and InDels.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @import plyr
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
incorporateGenomicVariants <- function(ProteoDiscography, aggregateSamples = FALSE, aggregateWithinExon = TRUE, aggregateWithinTranscript = TRUE, ignoreOverlappingMutations = TRUE, threads = 1) {

  # Input validation. -------------------------------------------------------

  checkmate::assertClass(ProteoDiscography, classes = "ProteoDiscography")
  checkmate::assertLogical(aggregateSamples)
  checkmate::assertLogical(aggregateWithinExon)
  checkmate::assertLogical(aggregateWithinTranscript)
  checkmate::assertLogical(ignoreOverlappingMutations)
  checkmate::assertNumeric(threads)

  # Check if any mutations have been imported.
  if(base::sum(S4Vectors::elementNROWS(ProteoDiscography@input.genomicVariants)) == 0) stop('No imported SNV, MNV or InDel events to incorporate.')

  ParallelLogger::logInfo('ProteoDisco - Adding SNV/MNV and InDel mutations to transcript sequences.')

  # Retrieve overlapping elements per mutation ------------------------------

  ParallelLogger::logTrace('\tProteoDisco - Finding overlap of mutations with CDS within the TxDb.')

  # Retrieve overlapping coding elements.
  allMuts <- GenomicRanges::sort(base::unique(base::unlist(ProteoDiscography@input.genomicVariants)))
  allOverlappingElements <- GenomicFeatures::cdsByOverlaps(ProteoDiscography@TxDb, allMuts, minoverlap = 1, type = 'any', column = c('CDSID', 'TXID', 'GENEID'))
  if(any(!lengths(allOverlappingElements$GENEID))) allOverlappingElements[which(!lengths(allOverlappingElements$GENEID)), ]$GENEID <- "Unknown"

  # Convert relevant information to tibble for faster indexing.
  allOverlappingElements.tibble <- tibble::as_tibble(S4Vectors::mcols(allOverlappingElements))

  # Determine overlapping CDS per mutation.
  perMutOverlappingElements <- IRanges::findOverlaps(query = allMuts, subject = allOverlappingElements)

  # Aggregate per mutation.
  perMutAggregate <- tibble::as_tibble(perMutOverlappingElements) %>%
    dplyr::group_by(.data$queryHits) %>%
    dplyr::summarise(
      dataOverlap = list(allOverlappingElements.tibble[.data$subjectHits,]),
      mutIndex = base::unique(.data$queryHits),
      CDSID = list(unique(.data$dataOverlap[[1]]$CDSID)),
      TXID = list(unique(.data$dataOverlap[[1]]$TXID)),
      GENEID = list(unique(.data$dataOverlap[[1]]$GENEID))
    ) %>%
    dplyr::select(c('mutIndex', 'CDSID', 'TXID', 'GENEID')) %>%
    dplyr::ungroup()

  # Merge overlapping genetic elements back to input mutations and generate a tibble.
  allMuts.Overlap <- allMuts[perMutAggregate$mutIndex,]
  S4Vectors::mcols(allMuts.Overlap) <- base::cbind(S4Vectors::mcols(allMuts.Overlap), perMutAggregate)

  allMuts <- tibble::as_tibble(allMuts.Overlap)
  allMuts <- allMuts %>% dplyr::select(c('seqnames', 'start', 'end', 'width',
                                         'strand', 'ref', 'alt',
                                         'mutationalType', 'sample', 'mutIndex',
                                         'CDSID', 'TXID', 'GENEID')) %>%
    dplyr::distinct()


  # Incorporate mutations into exons ----------------------------------------

  ParallelLogger::logInfo('\tProteoDisco - Incorporating mutations within the CDS sequence(s).')

  # If set to sample aggregation modus, simply set the name of all samples the same and only a single 'split' will be performed.
  if(aggregateSamples) allMuts$sample <- 'Aggregated'

  # Perform per sample (or on aggregate).
  mutantCDS.PerSample <- base::lapply(S4Vectors::split(allMuts, allMuts$sample), function(allMuts.PerSample){

    ParallelLogger::logTrace(sprintf('\t\tWorking on sample: %s (%s overlapping CDS; %s threads)', unique(allMuts.PerSample$sample), base::nrow(allMuts.PerSample), threads))

    # Generate tibble containing unique mutation:exon pairs.
    allMuts.PerSample$CDSID <- base::vapply(allMuts.PerSample$CDSID, paste0, collapse = ',', FUN.VALUE = c(''))
    allMuts.separatedOnCDS <- allMuts.PerSample %>% tidyr::separate_rows(.data$CDSID)

    # Get all relevant CDS.
    allCDS <- GenomicFeatures::cds(ProteoDiscography@TxDb, filter = list(cds_id = unique(allMuts.separatedOnCDS$CDSID)))

    if(length(allCDS) == 0) {
      ParallelLogger::logWarn(sprintf('\t\tNo CDS could be extracted for sample: %s', unique(allMuts.PerSample$sample)))
      return(character(0))
    }

    # Loop per exon harboring one or multiple overlapping mutations.
    # This returns an tibble with CDSID:mutantSequence which will be replacing the WT sequence of that exon in downstream functions.
    # Coding sequences are returned with correct orientation as they are transcribed into pre-mRNA.
    mutantCDS <- base::unlist(BiocParallel::bplapply(S4Vectors::split(allMuts.separatedOnCDS, allMuts.separatedOnCDS$CDSID), function(mutsPerCDS){

      ParallelLogger::logTrace(sprintf('\t\t\tWorking on CDS: %s (%s)', base::unique(mutsPerCDS$CDSID), base::unique(mutsPerCDS$sample)))

      # Add rownames as column to use for filtering purposes.
      mutsPerCDS <- mutsPerCDS %>% tibble::rownames_to_column()

      # Determine overlapping mutations and give appropriate response.
      if(aggregateWithinExon){
        mutsPerCDS.GRanges <- GenomicRanges::makeGRangesFromDataFrame(mutsPerCDS)

        # Reduce overlapping ranges.
        mutsPerCDS.GRanges.Reduced <- IRanges::reduce(mutsPerCDS.GRanges)

        # Per reduced ranges, determine overlap and only retain the first hit per reduced range by removing the others.
        overlappingMutsWithinCDS <- IRanges::findOverlaps(query = mutsPerCDS.GRanges, subject = mutsPerCDS.GRanges.Reduced, minoverlap = 1)

        # Remove self-hits which only overlap with themselves, i.e. 1 overlap count.
        selfHits <- IRanges::countOverlaps(mutsPerCDS.GRanges.Reduced, mutsPerCDS.GRanges, minoverlap = 1)
        selfHits <- which(selfHits == 1)
        if(length(selfHits) > 0) overlappingMutsWithinCDS <- overlappingMutsWithinCDS[S4Vectors::subjectHits(overlappingMutsWithinCDS) != selfHits,]

        # Determine number of overlapping mutations.
        overlappingMutsWithinCDSCount <- length(base::unique(S4Vectors::queryHits(overlappingMutsWithinCDS)))

        if(overlappingMutsWithinCDSCount > 0){
          if(ignoreOverlappingMutations){
            ParallelLogger::logWarn(sprintf('CDSID #%s has %s position-overlapping genomic mutations within the exon, only incorporating the first of these mutations based on the input order.', unique(mutsPerCDS$CDSID), overlappingMutsWithinCDSCount))

            # Remove non-first overlapping mutations.
            nonFirstHits <- base::sort(base::unique(base::unlist(base::lapply(S4Vectors::split(S4Vectors::queryHits(overlappingMutsWithinCDS), S4Vectors::subjectHits(overlappingMutsWithinCDS)), function(x) sort(x)[-1]))))
            mutsPerCDS <- mutsPerCDS %>% dplyr::filter(!.data$rowname %in% nonFirstHits)

          }else{
            ParallelLogger::logError(sprintf('CDSID #%s has %s position-overlapping genomic mutations within the exon and ignoreOverlappingMutations is set to FALSE. Ignoring this CDS', unique(mutsPerCDS$CDSID), overlappingMutsWithinCDSCount))
            return(NULL)
          }
        }
      }

      # Retrieve coding ranges within the supplied genome. We use CDS in order to have the 5' and 3' UTR removed.
      CDS.Granges <- allCDS[allCDS$cds_id == unique(mutsPerCDS$CDSID)]

      if(base::length(CDS.Granges) > 0){

        # Retrieve exonic sequence from supplied genomic sequences.
        CDS.SequenceWT <- BSgenome::getSeq(ProteoDiscography@genomeSeqs, CDS.Granges)[[1]]

        # Determine mutational positions within the exon, as if the start of the exon is position 1 (strand-specific).
        if(as.character(GenomicRanges::strand(CDS.Granges)) == '+'){
          mutsPerCDS$relativeCDSStart <- mutsPerCDS$start - IRanges::start(CDS.Granges) + 1
          mutsPerCDS$relativeCDSEnd <- mutsPerCDS$relativeCDSStart + (mutsPerCDS$width - 1)
          mutsPerCDS$refCDS <- mutsPerCDS$ref
          mutsPerCDS$refAlt <- mutsPerCDS$alt
        }else{
          mutsPerCDS$relativeCDSStart <- IRanges::end(CDS.Granges) - mutsPerCDS$end + 1
          mutsPerCDS$relativeCDSEnd <- mutsPerCDS$relativeCDSStart + mutsPerCDS$width - 1
          mutsPerCDS$refCDS <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mutsPerCDS$ref)))
          mutsPerCDS$refAlt <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mutsPerCDS$alt)))
        }

        # Determine which InDels delete the entire exon.
        if(base::any(mutsPerCDS$relativeCDSEnd > base::length(CDS.SequenceWT) & mutsPerCDS$relativeCDSStart < 1)){
          id <- mutsPerCDS[mutsPerCDS$relativeCDSEnd > base::length(CDS.SequenceWT) & mutsPerCDS$relativeCDSStart < 1,]
          mutsPerCDS$refCDS <- base::as.character(CDS.SequenceWT)
          mutsPerCDS$relativeCDSStart <- 1
          mutsPerCDS$relativeCDSEnd <- base::nchar(mutsPerCDS$refCDS)
          mutsPerCDS$refAlt <- ''
        }

        # Make boundary for exon-intron spanning mutations.
        if(base::any(mutsPerCDS$relativeCDSEnd > base::length(CDS.SequenceWT))){
          id <- mutsPerCDS$relativeCDSEnd > base::length(CDS.SequenceWT)
          mutsPerCDS[id,]$refCDS <- XVector::subseq(mutsPerCDS[id,]$refCDS, start = 1, end = mutsPerCDS[id,]$width - (mutsPerCDS[id,]$relativeCDSEnd - base::length(CDS.SequenceWT)))
          mutsPerCDS[id,]$relativeCDSEnd <- base::length(CDS.SequenceWT)
        }

        # Make boundary for intron-exon spanning mutations.
        if(base::any(mutsPerCDS$relativeCDSStart < 1)){
          id <- mutsPerCDS$relativeCDSStart < 1
          mutsPerCDS[id,]$refCDS <- XVector::subseq(mutsPerCDS[id,]$refCDS, start = abs(mutsPerCDS[id,]$relativeCDSStart) + 2, end = base::nchar(mutsPerCDS[id,]$refCDS))
          mutsPerCDS[id,]$refAlt <- XVector::subseq(mutsPerCDS[id,]$refAlt, start = abs(mutsPerCDS[id,]$relativeCDSStart) + 2, end = base::nchar(mutsPerCDS[id,]$refAlt))
          mutsPerCDS[id,]$relativeCDSStart <- 1
        }

        # Check if reference bases correspond.
        refCheck <- base::apply(mutsPerCDS, 1, function(mut){
          Biostrings::substr(CDS.SequenceWT, mut$relativeCDSStart, mut$relativeCDSEnd) == Biostrings::DNAString(mut$refCDS)
        })

        if(!all(refCheck)){
          stop(sprintf('Reference bases of %s mutations do not match reference bases of the overlapping exonic sequence (#%s) in the genome.', length(!refCheck), unique(mutsPerCDS$CDSID)))
        }

        # Replace the WT sequence with the mutation(s) on their relative positions.
        if(aggregateWithinExon){
          CDS.SequenceMut <- Biostrings::DNAStringSet(Biostrings::replaceAt(CDS.SequenceWT, at = IRanges::IRanges(mutsPerCDS$relativeCDSStart, mutsPerCDS$relativeCDSEnd), value = mutsPerCDS$refAlt))

          mutantSeq <- S4Vectors::DataFrame(
            CDS.SequenceMut = CDS.SequenceMut,
            CDS.ID = unique(mutsPerCDS$CDSID),
            GENE.ID = base::unique(base::unlist(mutsPerCDS$GENEID)),
            Tx.Strand = base::unique(base::as.character(GenomicRanges::strand(CDS.Granges))),
            mutations.genomic = base::paste(base::sprintf('g.%s:%s%s>%s', mutsPerCDS$seqnames, mutsPerCDS$start, mutsPerCDS$ref, mutsPerCDS$alt), collapse = ','),
            mutations.withinCDS = base::paste(base::sprintf('e.%s:%s%s>%s', unique(mutsPerCDS$CDSID), mutsPerCDS$relativeCDSStart, mutsPerCDS$refCDS, mutsPerCDS$refAlt), collapse = ',')
          )

          return(mutantSeq)

        } else {
          mutantSeq <- base::do.call(base::rbind, base::apply(mutsPerCDS, 1, function(mut){
            CDS.SequenceMut <- Biostrings::DNAStringSet(Biostrings::replaceAt(CDS.SequenceWT, at = IRanges::IRanges(mut$relativeCDSStart, mut$relativeCDSEnd), value = mut$refAlt))

            mutantSeq <- S4Vectors::DataFrame(
              CDS.SequenceMut = CDS.SequenceMut,
              CDS.ID = unique(mutsPerCDS$CDSID),
              GENE.ID = base::unique(base::unlist(mutsPerCDS$GENEID)),
              Tx.Strand = base::unique(base::as.character(GenomicRanges::strand(CDS.Granges))),
              mutations.genomic = base::paste(base::sprintf('g.%s:%s%s>%s', mutsPerCDS$seqnames, mutsPerCDS$start, mutsPerCDS$ref, mutsPerCDS$alt), collapse = ','),
              mutations.withinCDS = base::paste(base::sprintf('e.%s:%s%s>%s', unique(mutsPerCDS$CDSID), mutsPerCDS$relativeCDSStart, mutsPerCDS$refCDS, mutsPerCDS$refAlt), collapse = ',')
            )

            return(mutantSeq)

          }))
        }
        return(mutantSeq)

      } else {
        return(NULL)
      }
    }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, stop.on.error = TRUE, progressbar = TRUE)))

    # Return mutant exons per sample (or aggregate).
    mutantCDS.Combined <- base::do.call(base::rbind, mutantCDS)
    return(mutantCDS.Combined)

  })

  # Remove any samples without returnable (mutant) CDS
  mutantCDS.PerSample <- purrr::compact(mutantCDS.PerSample)

  # Generate mutant transcripts ---------------------------------------------

  ParallelLogger::logInfo('\tProteoDisco - Generating mutant transcript sequence(s) by replacing wild-type to mutant CDS.')

  # Per sample (or aggregate), generate the mutant transcripts.
  mutantTx.PerSample <- base::lapply(base::names(mutantCDS.PerSample), function(currentSample){

    # Get the mutations and mutant CDS of the sample.
    mutantCDS.currentSample <- mutantCDS.PerSample[base::names(mutantCDS.PerSample) == currentSample][[1]]

    ParallelLogger::logInfo(sprintf('\t\tWorking on sample: %s (%s mutant CDS; %s threads)', currentSample, base::nrow(mutantCDS.currentSample), threads))

    # Retrieve the transcripts of the overlapping CDS.
    Tx <- GenomicFeatures::transcripts(ProteoDiscography@TxDb, filter = list(cds_id = unique(mutantCDS.currentSample$CDS.ID)), column = 'TXID')

    # Retrieve all (incl. non-overlapping) CDS of all transcript.
    # We have to retrieve all TX:CDS because GenomicFeatures does not support pre-filtering..
    Tx.CDS <- BiocGenerics::unlist(GenomicFeatures::cdsBy(ProteoDiscography@TxDb, by = 'tx'))
    Tx.CDS$TXID <- base::names(Tx.CDS)

    # Filter on relevant TX.
    Tx.CDS <- Tx.CDS[Tx.CDS$TXID %in% base::unique(Tx$TXID)]
    Tx.CDS <- S4Vectors::split(Tx.CDS, Tx.CDS$TXID)

    # Per Tx, loop and generate the full mutant sequences by replacing the overlapping CDS with the mutant CDS sequences.
    mutSeqTx.CurrentSample <- BiocParallel::bplapply(Tx.CDS, function(Tx.CDS.Current) {

      ParallelLogger::logTrace(sprintf('\t\tWorking on TX: %s (%s)', BiocGenerics::unique(Tx.CDS.Current$TXID), currentSample))

      # Add the CDS as names which the DNAStringset inherits.
      base::names(Tx.CDS.Current) <- Tx.CDS.Current$cds_id

      # Retrieve the genomic sequences of each CDS.
      Tx.CDS.Current.SequenceWT <- BSgenome::getSeq(ProteoDiscography@genomeSeqs, Tx.CDS.Current)

      # Retrieve the mutant CDS falling within this Tx.
      mutantCDS.currentSample.currentTx <- mutantCDS.currentSample[mutantCDS.currentSample$CDS.ID %in% base::unique(Tx.CDS.Current$cds_id),]

      # Per mutant CDS, replace the WT CDS with the mutant CDS.
      # If aggregateWithinTranscript is set, simultaneously replace all mutant CDS within each Tx.
      # When multiple mutant-sequences for the same exon is present (with aggregateWithinTranscript = TRUE and aggregateWithinExon = FALSE),
      # only select the first viable option.
      # Else, generate separate Tx per mutant CDS.
      if(aggregateWithinTranscript){

        # Copy the WT sequence to tmp. object so we can overwrite it with mutant CDS.
        Tx.CDS.Current.SequenceWT.tmp <- Tx.CDS.Current.SequenceWT

        # Mutant sequence of current CDS.
        mutCDS <- S4Vectors::unique(mutantCDS.currentSample.currentTx$CDS.SequenceMut)
        names(mutCDS) <- S4Vectors::unique(mutantCDS.currentSample.currentTx$CDS.ID)

        # If only a single CDS per exon can be replaced.
        if(aggregateWithinExon & base::any(!base::duplicated(mutantCDS.currentSample.currentTx$CDS.ID))){

          # Replace WT to Mutant Seq.
          index <- which(base::names(Tx.CDS.Current.SequenceWT.tmp) %in% mutantCDS.currentSample$CDS.ID)
          Tx.CDS.Current.SequenceWT.tmp[index] <- mutCDS[names(Tx.CDS.Current.SequenceWT.tmp[index])]

          # Concatenate the CDS together into the full-length RNA-transcript.
          mutTx <- Biostrings::DNAStringSet(Biostrings::xscat(base::unlist(Tx.CDS.Current.SequenceWT.tmp)))

          # Add the mutant transcript sequence to the mutational information.
          mutantTx.currentSample.currentTx <- S4Vectors::DataFrame(
            Tx.Strand = base::unique(mutantCDS.currentSample.currentTx$Tx.Strand),
            mutations.genomic = base::paste(mutantCDS.currentSample.currentTx$mutations.genomic, collapse = ','),
            mutations.withinCDS = base::paste(mutantCDS.currentSample.currentTx$mutations.withinCDS, collapse = ','),
            numberOfMutationsInTx = base::sum(S4Vectors::elementNROWS(base::strsplit(mutantCDS.currentSample.currentTx$mutations.withinCDS, ','))),
            numberOfCDS = base::length(Tx.CDS.Current),
            numberOfMutCDS = base::length(base::unique(mutantCDS.currentSample.currentTx$CDS.ID)),
            geneID = base::unique(mutantCDS.currentSample.currentTx$GENE.ID),
            Tx.SequenceMut = mutTx
          )

          return(mutantTx.currentSample.currentTx)

        }else{

          # Multiple instances of the same exon need to be replaced within the same transcript.
          # Therefore generate multiple transcripts for these multi-exon replacements.
          ParallelLogger::logTrace(base::sprintf('\t\tThis TX has multiple forms of the same mutant-exon, therefore generating these as individual mutant-TX: %s (%s)', BiocGenerics::unique(Tx.CDS.Current$TXID), currentSample))

          # Loop over each mutant CDS.
          mutantTx.currentSample.currentTx <- base::do.call(base::rbind, base::lapply(base::seq_len(base::nrow(mutantCDS.currentSample.currentTx)), function(row){

            mutantCDS.currentSample.currentTx.currentRow <- mutantCDS.currentSample.currentTx[row,]

            # Copy the WT sequence to tmp. object so we can overwrite it with mutant CDS.
            Tx.CDS.Current.SequenceWT.tmp <- Tx.CDS.Current.SequenceWT

            # Mutant sequence of current CDS.
            mutCDS <- mutantCDS.currentSample.currentTx.currentRow$CDS.SequenceMut

            # Replace WT to Mutant Seq.
            Tx.CDS.Current.SequenceWT.tmp[base::names(Tx.CDS.Current.SequenceWT.tmp) %in% mutantCDS.currentSample.currentTx.currentRow$CDS.ID] <- mutCDS

            # Concetenate the CDS together into the full-length RNA-transcript.
            mutTx <- Biostrings::DNAStringSet(Biostrings::xscat(base::unlist(Tx.CDS.Current.SequenceWT.tmp)))

            # Subset the mutations we are currently incorporating.
            mutations.withinCDS <- base::unlist(Biostrings::strsplit(mutantCDS.currentSample.currentTx.currentRow$mutations.withinCDS, ','))
            mutations.genomic <- base::unlist(Biostrings::strsplit(mutantCDS.currentSample.currentTx.currentRow$mutations.genomic, ','))
            mutations.withinCDSInDex <- base::which(base::grepl(base::unique(mutantCDS.currentSample.currentTx.currentRow$CDS.ID), mutations.withinCDS))

            # Add the mutant transcript sequence to the mutational information.
            mutantTx.currentSample.currentTx <- S4Vectors::DataFrame(
              Tx.Strand = base::unique(mutantCDS.currentSample.currentTx.currentRow$Tx.Strand),
              mutations.genomic = base::paste(mutations.genomic[mutations.withinCDSInDex], collapse = ','),
              mutations.withinCDS = base::paste(mutations.withinCDS[mutations.withinCDSInDex], collapse = ','),
              numberOfMutationsInTx = base::sum(S4Vectors::elementNROWS(base::strsplit(mutations.withinCDS[mutations.withinCDSInDex], ','))),
              numberOfCDS = base::length(Tx.CDS.Current),
              numberOfMutCDS = base::length(base::unique(mutantCDS.currentSample.currentTx.currentRow$CDS.ID)),
              geneID = base::unique(mutantCDS.currentSample.currentTx$GENE.ID),
              Tx.SequenceMut = mutTx
            )

            return(mutantTx.currentSample.currentTx)

          }))

          return(mutantTx.currentSample.currentTx)

        }
      } else {

        # Loop over each mutant CDS.
        mutantTx.currentSample.currentTx <- base::do.call(base::rbind, base::lapply(base::seq_len(base::nrow(mutantCDS.currentSample.currentTx)), function(row){

          mutantCDS.currentSample.currentTx.currentRow <- mutantCDS.currentSample.currentTx[row,]

          # Copy the WT sequence to tmp. object so we can overwrite it with mutant CDS.
          Tx.CDS.Current.SequenceWT.tmp <- Tx.CDS.Current.SequenceWT

          # Mutant sequence of current CDS.
          mutCDS <- mutantCDS.currentSample.currentTx.currentRow$CDS.SequenceMut

          # Replace WT to Mutant Seq.
          Tx.CDS.Current.SequenceWT.tmp[base::names(Tx.CDS.Current.SequenceWT.tmp) %in% mutantCDS.currentSample.currentTx.currentRow$CDS.ID] <- mutCDS

          # Concetenate the CDS together into the full-length RNA-transcript.
          mutTx <- Biostrings::DNAStringSet(Biostrings::xscat(base::unlist(Tx.CDS.Current.SequenceWT.tmp)))

          # Subset the mutations we are currently incorporating.
          mutations.withinCDS <- base::unlist(Biostrings::strsplit(mutantCDS.currentSample.currentTx.currentRow$mutations.withinCDS, ','))
          mutations.genomic <- base::unlist(Biostrings::strsplit(mutantCDS.currentSample.currentTx.currentRow$mutations.genomic, ','))
          mutations.withinCDSInDex <- base::which(base::grepl(base::unique(mutantCDS.currentSample.currentTx.currentRow$CDS.ID), mutations.withinCDS))

          # Add the mutant transcript sequence to the mutational information.
          mutantTx.currentSample.currentTx <- S4Vectors::DataFrame(
            Tx.Strand = base::unique(mutantCDS.currentSample.currentTx.currentRow$Tx.Strand),
            mutations.genomic = base::paste(mutations.genomic[mutations.withinCDSInDex], collapse = ','),
            mutations.withinCDS = base::paste(mutations.withinCDS[mutations.withinCDSInDex], collapse = ','),
            numberOfMutationsInTx = base::sum(S4Vectors::elementNROWS(base::strsplit(mutations.withinCDS[mutations.withinCDSInDex], ','))),
            numberOfCDS = base::length(Tx.CDS.Current),
            numberOfMutCDS = base::length(base::unique(mutantCDS.currentSample.currentTx.currentRow$CDS.ID)),
            geneID = base::unique(mutantCDS.currentSample.currentTx$GENE.ID),
            Tx.SequenceMut = mutTx
          )

          return(mutantTx.currentSample.currentTx)
        }))

      }

      # Add sample name.
      mutantTx.currentSample.currentTx$sample <- currentSample

      return(mutantTx.currentSample.currentTx)

    }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE))

    # Add the TX.ID.
    mutSeqTx.CurrentSample.Combined <- do.call(rbind, unlist(mutSeqTx.CurrentSample))
    mutSeqTx.CurrentSample.Combined$Tx.ID <- base::rep(base::names(mutSeqTx.CurrentSample), S4Vectors::elementNROWS(mutSeqTx.CurrentSample))
    mutSeqTx.CurrentSample.Combined$sample <- currentSample

    return(mutSeqTx.CurrentSample.Combined)

  })

  # Combine all samples.
  mutantTx.allSamples <- base::do.call(base::rbind, mutantTx.PerSample)


  # Add 3' UTR to transcripts --------------------------------------------

  ParallelLogger::logTrace('Retrieving 3\' UTR')

  # Retrieve all 3' UTRs.
  allUTRs <- GenomicFeatures::threeUTRsByTranscript(ProteoDiscography@TxDb, use.names = FALSE)

  # Subset overlapping Tx.
  allUTRs <- base::unlist(allUTRs[names(allUTRs) %in% unique(mutantTx.allSamples$Tx.ID)])
  allUTRs$Tx.ID <- base::names(allUTRs)
  allUTRs$strand <- BiocGenerics::strand(allUTRs)

  # Retrieve the sequence.
  allUTRs.Seq <- Biostrings::getSeq(ProteoDiscography@genomeSeqs, allUTRs)

  # Collapse exon-exon spanning UTRs based on orientation.
  allUTRs.Seq <- tibble::tibble(seq.UTR = base::as.character(allUTRs.Seq), Tx.Id = base::names(allUTRs.Seq), strand = base::as.character(BiocGenerics::strand(allUTRs))) %>%
    dplyr::group_by(.data$Tx.Id) %>%
    dplyr::summarize(
      seq.UTR = base::ifelse(base::unique(.data$strand) == '+', base::paste0(.data$seq.UTR, collapse = ''), paste0(base::rev(.data$seq.UTR), collapse = ''))
    ) %>%
    dplyr::ungroup()

  # For Tx without 3' UTR in the TxDb, find the next 99nt. after last position of CDS.
  Tx.WithoutUTR <- mutantTx.allSamples[!mutantTx.allSamples$Tx.ID %in% allUTRs.Seq$Tx.Id,]

  if(base::nrow(Tx.WithoutUTR) > 0){

    Tx.WithoutUTR <- GenomicFeatures::transcripts(ProteoDiscography@TxDb, filter = list(tx_id = base::unique(Tx.WithoutUTR$Tx.ID)))

    # Add UTR to mut. TX-sequences (if UTR was present).
    Tx.WithoutUTR.Gr <- GenomicRanges::GRanges(
      seqnames = base::as.character(GenomeInfoDb::seqnames(Tx.WithoutUTR)),
      ranges = IRanges::IRanges(
        start = ifelse(BiocGenerics::strand(Tx.WithoutUTR) == '+', BiocGenerics::end(Tx.WithoutUTR) + 1, BiocGenerics::start(Tx.WithoutUTR) - 99),
        end = ifelse(BiocGenerics::strand(Tx.WithoutUTR) == '+', BiocGenerics::end(Tx.WithoutUTR) + 99, BiocGenerics::start(Tx.WithoutUTR) - 1)
      ),
      strand = BiocGenerics::strand(Tx.WithoutUTR)
    )

    base::names(Tx.WithoutUTR.Gr) <- Tx.WithoutUTR$tx_id

    # Retrieve 99nt 3' UTR.
    withoutUTRs.Seq <- Biostrings::getSeq(ProteoDiscography@genomeSeqs, Tx.WithoutUTR.Gr)
    withoutUTRs.Seq <- tibble::tibble(seq.UTR = as.character(withoutUTRs.Seq), Tx.Id = base::names(withoutUTRs.Seq))

    # Combine known and generated 3' UTRs.
    allUTRs.Seq <- dplyr::bind_rows(allUTRs.Seq, withoutUTRs.Seq)

  }

  # Match to input DataFrame.
  allUTRs.Seq <- allUTRs.Seq[base::match(mutantTx.allSamples$Tx.ID, allUTRs.Seq$Tx.Id),]

  # Combine TX and UTR.
  mutantTx.allSamples$Tx.SequenceMut <- Biostrings::xscat(mutantTx.allSamples$Tx.SequenceMut, allUTRs.Seq$seq.UTR)


  # Add variant transcripts to ProteoDiscography -----------------------------
  ParallelLogger::logInfo('\tProteoDisco - Adding mutant transcripts to the ProteoDiscography.')

  ProteoDiscography <- setMutantTranscripts(ProteoDiscography,
                                            transcripts = mutantTx.allSamples,
                                            slotType = 'genomicVariants')

  # Return statement --------------------------------------------------------

  return(ProteoDiscography)

}
