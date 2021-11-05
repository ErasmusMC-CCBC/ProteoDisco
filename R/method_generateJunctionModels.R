#' @title Generate putative transcript-models derived from splice-junctions.
#' @description Generates splicing-isoforms using the supplied splice-junctions (SJ) within the ProteoDiscography.
#'
#' Supplied junctions will be transformed into novel gene-models based on the nearest (or overlapping) exon with the TxDb
#' of the given {ProteoDiscography}; strand-information is taken into account. If no strand-information is given, the nearest (or overlapping) known
#' exon is assigned.
#'
#' It will derive the putative gene-model based on the assigned exons as annotated within the TxDb.
#' E.g., if matched with the second exon of geneA (+) and the fourth exon of geneB (+) it will generate the following gene-model:
#'
#' geneA-Exon1, geneA-Exon2, geneB-Exon4, geneB-Exon5, ...
#'
#' Users can also specify the max. search-window (in bp) in which the nearest canonical exonic boundary should fall.
#'
#' @param ProteoDiscography (\link[ProteoDisco]{ProteoDiscography}): ProteoDiscography object which stores the annotation and genomic sequences.
#' @param maxDistance (integer): Max. distance (>=0 bp) from splice-junction to nearest (or overlapping) exon in bp.
#' Setting this to high numbers will (erroneously) assign a distant exon as 'neighboring' resulting in unlikely models.
#' If the SJ is beyond this distance, it will be assigned as an cryptic exon and the length of this exon is set with the [maxCrypticSize] parameter.
#' @param maxCrypticSize (integer): The max. extension of bp (respective to orientation) beyond the SJ if it is assigned as a cryptic exon.
#' I.e., the nr. of bp that will be used in determining the putative transcript sequence for each cryptic junction.
#' @param skipCanonical (logical): Should canonical exon-exon junctions be skipped (TRUE) or generated (FALSE)?
#' @param threads (integer): Number of threads.
#'
#' @examples
#'
#'  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#'    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'  )
#'
#'  # Import splice-junctions (even spanning different chromosomes) based on our format.
#'  testSJ <- readr::read_tsv(system.file('extdata', 'validationSetSJ_hg19.txt', package = 'ProteoDisco'))
#'
#'  # Add custom SJ to ProteoDiscography.
#'  ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
#'    ProteoDiscography = ProteoDiscography.hg19,
#'    inputSpliceJunctions = testSJ
#'  )
#'
#'  # Generate junction-models from non-canonical splice-junctions.
#'  ProteoDiscography.hg19 <- ProteoDisco::generateJunctionModels(
#'    ProteoDiscography = ProteoDiscography.hg19,
#'    # Max. distance from a known exon-boundary before introducing a novel exon.
#'    # If an adjacent exon is found within this distance, it will shorten or elongate that exon towards the SJ.
#'    maxDistance = 150,
#'    # Should we skip known exon-exon junctions (in which both the acceptor and donor are located on known adjacent exons within the same transcript)
#'    skipCanonical = TRUE,
#'    # Perform on multiple threads (optional)
#'    threads = 1
#'  )
#'
#' @return {ProteoDiscography} with derived splice-isoforms.
#' @concept methods
#' @keywords methods
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @author Wesley van de Geer \email{w.vandegeer@erasmusmc.nl}
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
generateJunctionModels <- function(ProteoDiscography, maxDistance = 150, maxCrypticSize = 75, skipCanonical = TRUE, threads = 1){

  # Input validation. -------------------------------------------------------

  checkmate::assertClass(ProteoDiscography, classes = 'ProteoDiscography')
  checkmate::assertNumber(maxDistance, lower = 0, null.ok = FALSE)
  checkmate::assertNumber(maxCrypticSize, lower = 0, null.ok = FALSE)
  checkmate::assertLogical(skipCanonical)
  checkmate::assertNumeric(threads)

  sprintf('ProteoDisco - Generating splice-junctions models based on %s unique (sample-specific or aggregated) splice-junctions.', base::nrow(ProteoDiscography@input.spliceJunctions)) %>% ParallelLogger::logInfo()


  # Internal functions ------------------------------------------------------

  # Internal function to detect nearest / overlapping CDS.
  .findCDS <- function(TxDb.CDS, junctions, isJunctionA, maxDistance){

    sprintf('Retrieving %s-adjacent or overlapping CDS.', ifelse(isJunctionA, '5\'', '3\'')) %>% ParallelLogger::logTrace()

    # Convert splice-junction to GRanges.
    junction.GR <- data.frame(site = junctions) %>%
      tidyr::separate(col = 'site', into = c('chr', 'start', 'strand'), sep = ':', remove = FALSE) %>%
      dplyr::mutate(
        end = .data$start,
        cdsType = NA,
        cds_id = NA,
        chrDerivedCDS = NA,
        startAssignedCDS = NA,
        endAssignedCDS = NA,
        startDerivedCDS = NA,
        endDerivedCDS = NA,
        strandDerivedCDS = NA,
        distanceAssignedCDS = NA,
        junctionType = 'No known adjacent CDS.'
      ) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    # Clean-up nomenclature chromosomes.
    GenomeInfoDb::seqlevelsStyle(junction.GR) <- GenomeInfoDb::seqlevelsStyle(TxDb.CDS)

    # Detect if junction overlaps with an existing CDS.
    junction.GR$cds_id <- IRanges::findOverlaps(junction.GR, TxDb.CDS, select = 'first', type = 'within', maxgap = 0, ignore.strand = FALSE)
    if(any(!is.na(junction.GR$cds_id))) junction.GR[!is.na(junction.GR$cds_id),]$cdsType <- 'Overlap'

    # If no overlapping CDS is found, retrieve the nearest CDS.
    # Retrieve all adjacent CDS based on 5' -> 3' coordinates (i.e., in reverse for antisense genes).
    if(isJunctionA){
      if(any(GenomicRanges::strand(junction.GR) != '-')) junction.GR[GenomicRanges::strand(junction.GR) != '-' & is.na(junction.GR$cds_id)]$cds_id <- IRanges::follow(junction.GR[GenomicRanges::strand(junction.GR) != '-' & is.na(junction.GR$cds_id)], TxDb.CDS, select = 'last')
      if(any(GenomicRanges::strand(junction.GR) == '-')) junction.GR[GenomicRanges::strand(junction.GR) == '-' & is.na(junction.GR$cds_id)]$cds_id <- IRanges::precede(junction.GR[GenomicRanges::strand(junction.GR) == '-' & is.na(junction.GR$cds_id)], TxDb.CDS, select = 'first')

    }else{
      if(any(GenomicRanges::strand(junction.GR) != '-')) junction.GR[GenomicRanges::strand(junction.GR) != '-' & is.na(junction.GR$cds_id)]$cds_id <- IRanges::precede(junction.GR[GenomicRanges::strand(junction.GR) != '-' & is.na(junction.GR$cds_id)], TxDb.CDS, select = 'first')
      if(any(GenomicRanges::strand(junction.GR) == '-')) junction.GR[GenomicRanges::strand(junction.GR) == '-' & is.na(junction.GR$cds_id)]$cds_id <- IRanges::follow(junction.GR[GenomicRanges::strand(junction.GR) == '-' & is.na(junction.GR$cds_id)], TxDb.CDS, select = 'last')
    }

    # Annotate type of overlap.
    junction.GR[is.na(junction.GR$cdsType) & !is.na(junction.GR$cds_id),]$cdsType <- 'Nearest'
    if(any(is.na(junction.GR$cds_id))) junction.GR[is.na(junction.GR$cds_id),]$cdsType <- 'None'

    # Retrieve annotation of assigned CDS.
    junction.GR[!is.na(junction.GR$cds_id)]$chrDerivedCDS <- GenomeInfoDb::seqnames(junction.GR[!is.na(junction.GR$cds_id)])
    junction.GR[!is.na(junction.GR$cds_id)]$startAssignedCDS <- GenomicRanges::start(TxDb.CDS[match(junction.GR[!is.na(junction.GR$cds_id)]$cds_id, TxDb.CDS$cds_id)])
    junction.GR[!is.na(junction.GR$cds_id)]$endAssignedCDS <- GenomicRanges::end(TxDb.CDS[match(junction.GR[!is.na(junction.GR$cds_id)]$cds_id, TxDb.CDS$cds_id)])
    junction.GR[!is.na(junction.GR$cds_id)]$strandDerivedCDS <- GenomicRanges::strand(TxDb.CDS[match(junction.GR[!is.na(junction.GR$cds_id)]$cds_id, TxDb.CDS$cds_id)])

    # Derive spliced-CDS and calculate distance from start (if 5' junction) or end (if 3' junction) of overlapping CDS to splice-junction.
    if(isJunctionA){

      # Retrieve coordinates of derived CDS.
      # For 5' CDS, this is the start of the adjacent or overlapping CDS until the 5' SJ (-1 bp).
      junction.GR[!is.na(junction.GR$cds_id)]$startDerivedCDS <- junction.GR[!is.na(junction.GR$cds_id)]$startAssignedCDS
      junction.GR[!is.na(junction.GR$cds_id)]$endDerivedCDS <- GenomicRanges::start(junction.GR[!is.na(junction.GR$cds_id)]) - 1
      junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS <- junction.GR[!is.na(junction.GR$cds_id)]$endAssignedCDS - junction.GR[!is.na(junction.GR$cds_id)]$endDerivedCDS

      # Annotate type of junction (If directly flanking, it's a known junction).
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS == 0, 'Known 5\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS > 0, 'Shortening of known 5\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS < 0, 'Extension of known 5\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)

    }else{

      # Retrieve coordinates of derived CDS.
      # For 3' CDS, this is the SJ (+1bp) until the end of the 3'-adjacent CDS.
      junction.GR[!is.na(junction.GR$cds_id)]$startDerivedCDS <- GenomicRanges::start(junction.GR[!is.na(junction.GR$cds_id)]) + 1
      junction.GR[!is.na(junction.GR$cds_id)]$endDerivedCDS <- GenomicRanges::end(TxDb.CDS[match(junction.GR[!is.na(junction.GR$cds_id)]$cds_id, TxDb.CDS$cds_id)])
      junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS <- junction.GR[!is.na(junction.GR$cds_id)]$startAssignedCDS - junction.GR[!is.na(junction.GR$cds_id)]$startDerivedCDS

      # Annotate type of junction (If directly flanking, it's a known junction).
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS == 0, 'Known 3\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS > 0, 'Extension of known 3\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)
      junction.GR[!is.na(junction.GR$cds_id)]$junctionType <- ifelse(junction.GR[!is.na(junction.GR$cds_id)]$distanceAssignedCDS < 0, 'Shortening of known 3\' CDS-junction', junction.GR[!is.na(junction.GR$cds_id)]$junctionType)

    }

    # Set SJ with their distance > maxDistance as cryptic CDSs rather than extension / shortenings of existing CDS.
    idNonKnownSJ <- grepl('Extension|Shortening', junction.GR$junctionType) & (abs(junction.GR$distanceAssignedCDS) > maxDistance)
    if(sum(idNonKnownSJ) > 0) junction.GR[idNonKnownSJ,]$junctionType <- sprintf('Possible %s cryptic CDS due to distance from known junction', ifelse(isJunctionA, '5\'', '3\''))

    # Convert to tibble.
    junction.DF <- tibble::as_tibble(GenomicRanges::mcols(junction.GR))

    # Add prefix to colnames.
    if(isJunctionA) colnames(junction.DF) <- paste0('JunctionA.', colnames(junction.DF))
    if(!isJunctionA) colnames(junction.DF) <- paste0('JunctionB.', colnames(junction.DF))

    # Convert to tibble and return.
    return(junction.DF)

  }


  # Retrieve CDS from TxDb --------------------------------------------------

  ParallelLogger::logTrace('Retrieving all CDS-records from the TxDb.')
  TxDb.CDS <- GenomicFeatures::cds(ProteoDiscography@TxDb, columns = c('cds_id', 'gene_id', 'tx_id', 'exon_rank'))
  TxDb.CDS <- GenomicRanges::sort.GenomicRanges(TxDb.CDS)


  # Retrieve nearest / overlapping 5' and 3' CDS ----------------------------

  # Retrieve imported splice-junctions.
  data.SpliceJunctions <- ProteoDiscography@input.spliceJunctions

  data.SpliceJunctions.perSample <- base::split(data.SpliceJunctions, data.SpliceJunctions$sample)

  # Loop per sample.
  data.SpliceJunctions.perSample <- BiocParallel::bplapply(data.SpliceJunctions.perSample, function(x){

    ParallelLogger::logTrace('Retrieving overlapping or adjacent CDS for each splice-junction: ', unique(x$sample))

    # Retrieve adjacent or overlapping 5' and 3' CDSs.
    junctionA <- .findCDS(TxDb.CDS, junctions = x$junctionA, isJunctionA = TRUE, maxDistance = maxDistance)
    junctionB <- .findCDS(TxDb.CDS, junctions = x$junctionB, isJunctionA = FALSE, maxDistance = maxDistance)

    # Combine the information from both junctions.
    x <- x %>%
      dplyr::left_join(junctionA, by = c('junctionA' = 'JunctionA.site')) %>%
      dplyr::left_join(junctionB, by = c('junctionB' = 'JunctionB.site')) %>%
      dplyr::distinct()

    # Return.
    return(x)

  }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, stop.on.error = TRUE, progressbar = TRUE))


  # Generate potential splice-isoform ---------------------------------------

  # Loop per sample.
  data.SpliceTranscripts.perSample <- BiocParallel::bplapply(data.SpliceJunctions.perSample, function(sample.SJ){

    ParallelLogger::logTrace(sprintf('Generating splice-transcripts from %s SJ: %s', base::nrow(sample.SJ), base::unique(sample.SJ$sample)))

    # Remove SJ prior to any known CDS or after any known CDS.
    sample.SJ <- sample.SJ %>% dplyr::filter(.data$JunctionA.cdsType != 'None' & .data$JunctionB.cdsType != 'None')

    ## Loop per SJ ----

    sample.SJ.withSeq <- dplyr::bind_rows(BiocParallel::bplapply(base::seq_len(base::nrow(sample.SJ)), function(i){

      SJ <- sample.SJ[i,]

      ### Retrieve the upstream and downstream CDS sequences ----

      cds.A <- TxDb.CDS[TxDb.CDS$cds_id == as.integer(SJ$JunctionA.cds_id)]
      cds.B <- TxDb.CDS[TxDb.CDS$cds_id == as.integer(SJ$JunctionB.cds_id)]

      ### Check if SJ forms a pre-existing canonical junction. ----

      if(skipCanonical){

        # Check if Tx match and the exon-ranks are following.
        if(base::any(base::unlist(DelayedArray::unique(c(BiocGenerics::paste(cds.A$tx_id, cds.A$exon_rank + 1), BiocGenerics::paste(cds.A$tx_id, cds.A$exon_rank - 1)))) %in% base::unlist(BiocGenerics::paste(cds.B$tx_id, cds.B$exon_rank)))){

          ParallelLogger::logTrace('Skipping due to canonical status: ', base::ifelse(!is.na(SJ$identifier), SJ$identifier, sprintf('A:%s;B:%s', SJ$junctionA, SJ$junctionB)))

          SJ$identifierFASTA <- NA
          SJ$Tx.SequenceMut <- NA
          SJ$AA.SequenceMut <- NA

          # Return empty sequences.
          return(SJ)

        }
      }

      # Define CDS A.
      if(!base::grepl('cryptic', SJ$JunctionA.junctionType)){

        # Replace non-cryptic CDS with derived coordinates from near-adjacent CDS.
        IRanges::ranges(cds.A) <- IRanges::IRanges(
          start = SJ$JunctionA.startDerivedCDS,
          end = SJ$JunctionA.endDerivedCDS
        )

      }else{

        # Check if SJ is not > maxDistance from assigned CDS.
        distanceSJtoCDS <- base::abs((as.integer(Biostrings::strsplit(SJ$junctionA, ':')[[1]][2]) - 1) - IRanges::start(cds.A))

        if(distanceSJtoCDS <= maxDistance){

          # Replace cryptic CDS with SJ-coordinates + user-given maxCrypticSize.
          IRanges::ranges(cds.A) <- IRanges::IRanges(
            start = (as.integer(Biostrings::strsplit(SJ$junctionA, ':')[[1]][2]) - 1) - maxCrypticSize,
            end = as.integer(Biostrings::strsplit(SJ$junctionA, ':')[[1]][2]) - 1
          )

        }else{
          sprintf('Skipping: SJ > maxDistance from assigned CDS: %s (%s - %s)', SJ$identifier, SJ$junctionA, SJ$junctionB) %>% ParallelLogger::logTrace()

          SJ$identifierFASTA <- NA
          SJ$Tx.SequenceMut <- NA
          SJ$AA.SequenceMut <- NA

          # Return empty sequences.
          return(SJ)
        }
      }

      # Define CDS B.
      if(!base::grepl('cryptic', SJ$JunctionB.junctionType)){

        # Replace non-cryptic CDS with derived coordinates from near-adjacent CDS.
        IRanges::ranges(cds.B) <- IRanges::IRanges(
          start = SJ$JunctionB.startDerivedCDS,
          end = SJ$JunctionB.endDerivedCDS
        )
      }else{

        # Check if SJ is not > maxDistance from assigned CDS.
        distanceSJtoCDS <- base::abs(as.integer(Biostrings::strsplit(SJ$junctionB, ':')[[1]][2]) + 1 - IRanges::start(cds.B))

        if(distanceSJtoCDS <= maxDistance){

          # Replace cryptic CDS with SJ-coordinates + user-given maxCrypticSize.
          IRanges::ranges(cds.B) <- IRanges::IRanges(
            start = as.integer(Biostrings::strsplit(SJ$junctionB, ':')[[1]][2]) + 1,
            end = GenomicRanges::start(cds.B) + maxCrypticSize
          )
        }else{
          sprintf('Skipping: SJ > maxDistance from assigned CDS: %s (%s - %s)', SJ$identifier, SJ$junctionA, SJ$junctionB) %>% ParallelLogger::logTrace()

          SJ$identifierFASTA <- NA
          SJ$Tx.SequenceMut <- NA
          SJ$AA.SequenceMut <- NA

          # Return empty sequences.
          return(SJ)
        }

      }

      # Retrieve the sequences of the two CDS (already takes care of strand).
      cdsSeq.A <- Biostrings::getSeq(ProteoDiscography@genomeSeqs, cds.A)
      cdsSeq.B <- Biostrings::getSeq(ProteoDiscography@genomeSeqs, cds.B)

      # Perform three-frame translation for +/+ and -/- junctions (based on strand).
      if(all(GenomicRanges::strand(cds.A) == GenomicRanges::strand(cds.B))){

        # Combine sequences (A -> B)
        if(all(GenomicRanges::strand(cds.A) == '+')) cdsSeq.Combined <- Biostrings::xscat(cdsSeq.A, cdsSeq.B)

        # Combine sequences (A <- B)
        if(all(GenomicRanges::strand(cds.A) == '-')) cdsSeq.Combined <- Biostrings::xscat(cdsSeq.B, cdsSeq.A)

        # Three frame-translation (A -> B (+/+) or B <- A (-/-))
        # Cut mutant sequence(s) to first terminal.
        cdsSeq.Combined.AA <- base::list(base::unlist(BiocGenerics::lapply(base::seq_len(3), function(p) base::gsub('\\*', '', as.character(base::suppressWarnings(Biostrings::translate(Biostrings::subseq(cdsSeq.Combined, p))))))))

        SJ$Tx.SequenceMut <- base::as.character(cdsSeq.Combined)
        SJ$AA.SequenceMut <- cdsSeq.Combined.AA

      }else{

        # Perform six-frame translation if +/- / -/+ junctions.
        cdsSeq.Combined <- Biostrings::xscat(cdsSeq.A, cdsSeq.B)

        # Six-frame translation.
        # Cut mutant sequence(s) to first terminal.
        cdsSeq.Combined.AA <- base::unlist(BiocGenerics::lapply(base::seq_len(3), function(p) base::gsub('\\*', '', base::as.character(base::suppressWarnings(Biostrings::translate(Biostrings::subseq(c(cdsSeq.Combined, Biostrings::reverseComplement(cdsSeq.Combined)), p)))))))

        SJ$Tx.SequenceMut <- base::as.character(cdsSeq.Combined)
        SJ$AA.SequenceMut <- base::list(cdsSeq.Combined.AA)

      }

      # Generate identifier.
      SJ$identifierFASTA <- sprintf(
        '%sA:%s[%s];B:%s[%s]|A(%s):%s;B(%s):%s',
        ifelse(!is.na(SJ$identifier), sprintf('%s|', SJ$identifier), ''),
        BiocGenerics::paste(unique(cds.A$gene_id), collapse = ', '),
        cds.A$cds_id,
        BiocGenerics::paste(unique(cds.B$gene_id), collapse = ', '),
        cds.B$cds_id,
        SJ$JunctionA.junctionType,
        BiocGenerics::paste(GenomicRanges::seqnames(cds.A), GenomicRanges::start(cds.A), GenomicRanges::end(cds.A), GenomicRanges::strand(cds.A), sep = ':'),
        SJ$JunctionB.junctionType,
        BiocGenerics::paste(GenomicRanges::seqnames(cds.B), GenomicRanges::start(cds.B), GenomicRanges::end(cds.B), GenomicRanges::strand(cds.B), sep = ':')

      )

      # Return the SJ with putative AA-sequences.
      return(SJ)

    }, BPPARAM = BiocParallel::MulticoreParam(workers = threads, stop.on.error = TRUE, progressbar = TRUE)))

    # Print message.
    base::sprintf('%s - Out of %s given SJ, we converted %s SJ to splice-transcripts%s.', unique(sample.SJ.withSeq$sample), base::nrow(sample.SJ.withSeq), base::nrow(sample.SJ.withSeq %>% dplyr::filter(!is.na(.data$Tx.SequenceMut))), ifelse(skipCanonical, ' (non-canonical only)', '')) %>% ParallelLogger::logInfo()

    # Filter out non-used SJ.
    sample.SJ.withSeq <- sample.SJ.withSeq %>% dplyr::filter(!is.na(.data$Tx.SequenceMut))

    # Convert to DataFrame.
    sample.SJ.withSeq <- S4Vectors::DataFrame(sample.SJ.withSeq)

    # Return statement.
    return(sample.SJ.withSeq)

  })

  # Combine all samples.
  data.SpliceTranscripts.combined <- BiocGenerics::do.call(BiocGenerics::rbind, data.SpliceTranscripts.perSample)

  # Add generated SJ-transcript to ProteoDisco ------------------------------

  ProteoDiscography <- ProteoDisco::setMutantTranscripts(x = ProteoDiscography, transcripts = data.SpliceTranscripts.combined, slotType = 'spliceJunctions')


  # Return statement --------------------------------------------------------

  return(ProteoDiscography)

}
