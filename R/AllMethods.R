
# Getter methods ----------------------------------------------------------

#' @rdname getDiscography
#' @exportMethod getDiscography
setMethod('getDiscography', 'ProteoDiscography', function(x){
  
  # Retrieve imported records -----------------------------------------------
  
  y <- list(
    genomicVariants = x@input.genomicVariants,
    spliceJunctions = x@input.spliceJunctions,
    manualSequences = x@input.manualSequences
  )
  
  # Return list.
  return(y)
  
})

#' @rdname summary
#' @exportMethod summary
setMethod('summary', 'ProteoDiscography', function(object, verbose = TRUE){
  
  ### Overview of imported mutations -----------------------------------------
  
  # Total nr. of imported mutations per sample.
  if(base::length(object@input.genomicVariants) > 0){
    
    sample.ImportedMuts <- dplyr::bind_rows(base::lapply(object@input.genomicVariants, function(x){
      
      tibble::as_tibble(S4Vectors::mcols(x)) %>% 
        dplyr::group_by(.data$mutationalType) %>% 
        dplyr::summarise(totalN = dplyr::n(), sample = base::unique(sample)) %>% 
        dplyr::ungroup()
      
    }))
    
  }else{
    sample.ImportedMuts <- tibble::tibble()
  }
  
  # Total nr. of imported manual sequences.
  if(base::nrow(object@input.manualSequences) > 0){
    sample.ImportedManualSequences <- tibble::as_tibble(object@input.manualSequences) %>% dplyr::group_by(sample) %>% dplyr::summarise(totalN = dplyr::n(), mutationalType = 'Manual sequences') %>% dplyr::ungroup()  
  }else{
    sample.ImportedManualSequences <- tibble::tibble()
  }
  
  # Total nr. of predicted exon-exon events.
  if(base::nrow(object@input.spliceJunctions) > 0){
    sample.spliceJunctions <- tibble::as_tibble(object@input.spliceJunctions) %>% dplyr::group_by(sample) %>% dplyr::summarise(totalN = dplyr::n(), mutationalType = 'Splice-junctions') %>% dplyr::ungroup()  
  }else{
    sample.spliceJunctions <- tibble::tibble()
  }
  
  sample.ImportedOverview <- dplyr::bind_rows(sample.ImportedMuts, sample.ImportedManualSequences, sample.spliceJunctions)
  
  # Transpose.
  if(base::nrow(sample.ImportedOverview) > 0) sample.ImportedOverview <- sample.ImportedOverview %>% tidyr::spread(key = .data$mutationalType, value = .data$totalN)
  
  # Get the metadata of the TxDb.
  meta.txdb <- S4Vectors::metadata(object@TxDb)
  
  # Total number of unique samples.
  sampleTable <- 0
  if(base::nrow(sample.ImportedOverview) > 0) sampleTable <- dplyr::n_distinct(sample.ImportedOverview$sample)
  
  # Total number of distinct (mutated and manual) transcript sequences. 
  transTable <- base::sum(
    base::nrow(object@mutantTranscripts.genomicVariants),
    base::nrow(object@mutantTranscripts.spliceJunctions),
    base::nrow(object@input.manualSequences)
  )
  
  # Information on the Discography samples.
  if(verbose) cat(sprintf('This ProteoDiscography was initialized on %s for %s (%s).\n', 
                          object@metadata$CreatedOn, unique(organism(object)), 
                          unique(GenomeInfoDb::genome(object@genomeSeqs))))
  if(verbose) cat(sprintf('The underlying TxDb contains %s transcripts with %s exons.\n', 
                          meta.txdb[meta.txdb$name == 'transcript_nrow',]$value, 
                          meta.txdb[meta.txdb$name == 'exon_nrow',]$value))
  
  # Information on the input samples.
  if(base::nrow(sample.ImportedMuts) > 0){
    if(verbose) cat(sprintf('\nCurrently, %s genomic mutations from %s somatic (VCF or MAF) file(s) are imported.\n', base::sum(sample.ImportedMuts$totalN), dplyr::n_distinct(sample.ImportedMuts$sample)))
  } else {
    if(verbose) cat('\nNo VCF or MAF files containing SNV, InDels and/or MNV imported yet.\n')
  }
  
  # Information on prelim. overlap.
  if(!is.na(object@metadata$overlapUniqueCDS)){
    if(verbose) cat(sprintf('The imported genomic mutations overlap with %s distinct coding sequences from %s distinct genes.\n', object@metadata$overlapUniqueCDS, object@metadata$overlapUniqueGenes))
  }
  
  # Information of analyzed samples.
  if(!is.null(transTable) & nrow(sample.ImportedOverview) > 0) {
    if(verbose) cat(sprintf('\nMutational events (incl. manual sequences) from %s sample(s) have been processed; this generated %s mutant transcripts.\n', sampleTable, transTable))
  } else {
    if(verbose) cat('\nNo transcripts have been gathered.\n')
  }
  
  if(verbose) cat(sprintf('\nProteoDisco version during generation: v%s \n\n', object@metadata$ProteoDiscoVersion))
  
  return(list(overviewMutations = sample.ImportedOverview))
  
})

#' @rdname show
#' @exportMethod show
setMethod('show', 'ProteoDiscography', function(object){
  return(summary(object))
})

#' @rdname print
#' @exportMethod print
setMethod('print', 'ProteoDiscography', function(object){
  return(summary(object))
})

#' @rdname seqinfo
#' @exportMethod seqinfo
setMethod('seqinfo', 'ProteoDiscography', function(x){
  if(is.null(x@TxDb)) return('First add the TxDb to ProteoDiscography.')
  return(GenomicFeatures::seqinfo(x@TxDb))
})

#' @rdname seqlevels
#' @exportMethod seqlevels
setMethod('seqlevels', 'ProteoDiscography', function(x){
  if(is.null(x@TxDb)) return('First add the TxDb to ProteoDiscography.')
  return(GenomeInfoDb::seqlevels(x@TxDb))
})

#' @rdname organism
#' @exportMethod organism
setMethod('organism', 'ProteoDiscography', function(object){
  if(is.null(object@TxDb)) return('First add the TxDb to ProteoDiscography.')
  return(GenomicFeatures::organism(object@TxDb))
})

# Setter methods ----------------------------------------------------------


setMethod('.setGenomicVariants', 'ProteoDiscography', function(x, variants, removeExisting, overwriteDuplicateSamples){
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkClass(x, classes = 'ProteoDiscography')
  checkmate::checkClass(variants, classes = 'VRangesList')
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(overwriteDuplicateSamples)
  
  
  # Function ----------------------------------------------------------------
  
  # If set, remove all existing mutations within the ProteoDiscography.
  if(removeExisting) {
    x@input.genomicVariants <- VariantAnnotation::VRangesList()
  }
  
  # Overwrite duplicate samples or append novel mutations to existing samples.
  if(any(names(variants) %in% names(x@input.genomicVariants))){
    
    # Indexes of duplicate(s).
    duplicates <- which(names(variants) %in% names(x@input.genomicVariants))
    
    # Overwrite the old duplicates.
    if(overwriteDuplicateSamples){
      ParallelLogger::logTrace(sprintf('ProteoDisco - Found duplicate sample(s) - Overwriting: %s.', paste(unique(names(variants[duplicates])), collapse = ', ')))
      x@input.genomicVariants[duplicates] <- NULL
    }else{
      stop(sprintf('ProteoDisco - Found duplicate sample(s) (overwriteDuplicateSamples = FALSE): %s.', paste(unique(names(variants[duplicates])), collapse = ', ')))
    }
  }
  
  # Add new samples to VRangesList
  x@input.genomicVariants <- c(x@input.genomicVariants, variants)
  
  
  # Return ------------------------------------------------------------------
  
  return(x)
})

setMethod('.setManualSequences', 'ProteoDiscography', function(x, transcripts, removeExisting, overwriteDuplicateSamples) {
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkClass(x, classes = 'ProteoDiscography')
  checkmate::checkClass(transcripts, classes = 'DataFrame')
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(overwriteDuplicateSamples)
  
  
  # Function ----------------------------------------------------------------
  
  # If set, remove all existing mutations within the ProteoDiscography.
  if(removeExisting) {
    x@input.manualSequences <- S4Vectors::DataFrame()
  }
  
  # Overwrite duplicate samples or append sequences to existing samples.
  if(any(unique(transcripts$sample) %in% x@input.manualSequences$sample)){
    
    # Indexes of duplicate(s).
    allSamples <- unique(transcripts$sample)
    duplicates <- allSamples[which(allSamples %in% x@input.manualSequences$sample)]
    
    # Overwrite the old duplicates.
    if(overwriteDuplicateSamples){
      ParallelLogger::logTrace(sprintf('ProteoDisco - Found duplicate sample(s) - Overwriting: %s.', paste(duplicates, collapse = ', ')))
      x@input.manualSequences <- x@input.manualSequences[!x@input.manualSequences$sample %in% duplicates,]
    }else{
      stop(sprintf('ProteoDisco - Found duplicate sample(s) (overwriteDuplicateSamples = FALSE): %s.', paste(duplicates, collapse = ', ')))
    }
  }
  
  # Add new transcript sequences to DataFrame.
  transcripts <- base::rbind(x@input.manualSequences, transcripts)
  x@input.manualSequences <- transcripts
  
  # Add to ProteoDiscography.
  x <- setMutantTranscripts(x, transcripts, 'manualSequences')
  
  
  # Return ------------------------------------------------------------------
  
  return(x)
})


setMethod('.setSplicingJunctions', 'ProteoDiscography', function(x, spliceJunctions, removeExisting, overwriteDuplicateSamples) {
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkClass(x, classes = 'ProteoDiscography')
  checkmate::checkClass(spliceJunctions, classes = 'tbl_df')
  checkmate::assertLogical(removeExisting)
  checkmate::assertLogical(overwriteDuplicateSamples)
  
  
  # Function ----------------------------------------------------------------
  
  # If set, remove all existing mutations within the ProteoDiscography.
  if(removeExisting) {
    x@input.spliceJunctions <- tibble::tibble()
  }
  
  # Overwrite duplicate samples or append sequences to existing samples.
  if(nrow(x@input.spliceJunctions) > 0){
    if(any(unique(spliceJunctions$sample) %in% x@input.spliceJunctions$sample)){
      
      # Indexes of duplicate(s).
      allSamples <- unique(x@input.spliceJunctions$sample)
      duplicates <- allSamples[which(allSamples %in% x@input.spliceJunctions$sample)]
      
      # Overwrite the old duplicates.
      if(overwriteDuplicateSamples){
        ParallelLogger::logTrace(sprintf('ProteoDisco - Found duplicate sample(s) - Overwriting: %s.', paste(duplicates, collapse = ', ')))
        x@input.spliceJunctions <- x@input.spliceJunctions %>% dplyr::filter(!sample %in% duplicates)
      }else{
        stop(sprintf('ProteoDisco - Found duplicate sample(s) (overwriteDuplicateSamples = FALSE): %s.', paste(duplicates, collapse = ', ')))
      }
    }
  }
  
  # Add new transcript sequences to tibble
  spliceJunctions <- dplyr::bind_rows(x@input.spliceJunctions, spliceJunctions)
  spliceJunctions <- base::unique(spliceJunctions)
  
  # Add to ProteoDiscography.
  x@input.spliceJunctions <- spliceJunctions
  
  # Return ------------------------------------------------------------------
  
  return(x)
})


#' @exportMethod setMutantTranscripts
#' @rdname setMutantTranscripts
setMethod('setMutantTranscripts', 'ProteoDiscography', function(x, transcripts, slotType) {
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkClass(x, classes = 'ProteoDiscography')
  checkmate::checkClass(transcripts, classes = 'DataFrame')
  checkmate::checkCharacter(slotType, pattern = 'genomicVariants|manualSequences|spliceJunctions')
  
  if(!slotType %in% c('genomicVariants', 'manualSequences', 'spliceJunctions')){
    stop(sprintf('This ProteoDiscography slot has not yet been implemented: %s', slotType))
  }
  
  
  # Function ----------------------------------------------------------------
  
  translateSeq <- function(seq, geneticCode){
    
    if(methods::is(seq, 'character')) seq <- Biostrings::DNAStringSet(seq)
    
    base::suppressWarnings(
      AA.SequenceMut <- 
        # Translate
        seq %>% 
        Biostrings::translate(if.fuzzy.codon = 'solve', genetic.code = geneticCode) %>% 
        # Cutoff at first stop codon (*)
        base::gsub('\\*.*', '', .) %>% 
        # Convert to AAStringSet
        Biostrings::AAStringSet(.)
    )
    return(AA.SequenceMut)
  }
  
  if(slotType ==  'genomicVariants'){
    checkmate::checkClass(transcripts, classes = 'DataFrame')
    x@mutantTranscripts.genomicVariants <- transcripts
    
    # Translate to peptide sequences.
    x@mutantTranscripts.genomicVariants$AA.SequenceMut <- translateSeq(x@mutantTranscripts.genomicVariants$Tx.SequenceMut, x@GENETIC_CODE)
    
    return(x)
  }
  
  if(slotType ==  'manualSequences'){
    checkmate::checkClass(transcripts, classes = 'DataFrame')
    x@mutantTranscripts.manualSequences <- transcripts
    
    # Convert to peptide sequences.
    x@mutantTranscripts.manualSequences$AA.SequenceMut <- translateSeq(x@mutantTranscripts.manualSequences$Tx.SequenceMut, x@GENETIC_CODE)
    
    return(x)
  }
  
  if(slotType ==  'spliceJunctions'){
    checkmate::checkClass(transcripts, classes = 'DataFrame')
    x@mutantTranscripts.spliceJunctions <- transcripts
    
    return(x)
  }
  
})

#' @exportMethod setMutantTranscripts
#' @rdname setMutantTranscripts
setMethod('mutantTranscripts', 'ProteoDiscography', function(x) {
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkClass(x, classes = 'ProteoDiscography')
  
  # Function ----------------------------------------------------------------
  
  mutantTranscripts <- list()
  mutantTranscripts$genomicVariants <- x@mutantTranscripts.genomicVariants
  mutantTranscripts$spliceJunctions <- x@mutantTranscripts.spliceJunctions
  mutantTranscripts$manualSequences <- x@mutantTranscripts.manualSequences
  
  
  # Return statement --------------------------------------------------------
  
  return(mutantTranscripts)
  
})

#' @rdname checkProteotypicFragments
#' @exportMethod checkProteotypicFragments
setMethod('checkProteotypicFragments', 'ProteoDiscography', function(x, enzymUsed = 'Trypsin', missedCleavages = 0, additionalPeptides = NULL, checkWithinMutantSeqs = FALSE){
  
  # Input validation. -------------------------------------------------------
  
  checkmate::assertClass(x, classes = 'ProteoDiscography')
  checkmate::assertCharacter(enzymUsed)
  checkmate::assertInt(missedCleavages)
  checkmate::assertClass(additionalPeptides, classes = 'AAStringSet', null.ok = TRUE)
  checkmate::assertLogical(checkWithinMutantSeqs)
  
  # Check if any transcripts have been generated.
  if(base::sum(S4Vectors::elementNROWS(mutantTranscripts(x))) == 0) stop('No transcripts to check for proteotypic fragments! Please add / generate (mutant) transcripts first.')
  
  
  # Translate the TxDb. -----------------------------------------------------
  
  ParallelLogger::logInfo('ProteoDisco - Translating and cleaving the TxDb (and given peptide sequences) to generate proteotypic fragments. This can take a moment for large TxDb.')
  
  # Cleave the comparison database and only keep distinct records.
  TxDb.Fragments <- GenomicFeatures::cdsBy(x@TxDb, by = 'tx', use.names = TRUE) %>% 
    # Extract the DNA sequence based on CDS.
    GenomicFeatures::extractTranscriptSeqs(x = x@genomeSeqs) %>% 
    # Translate using the specified genetic code.
    Biostrings::translate(genetic.code = x@GENETIC_CODE, if.fuzzy.codon = 'solve') %>% 
    # Add the user-specific peptide-sequences.
    c(additionalPeptides) %>% 
    # Cleave using the specified settings.
    cleaver::cleave(enzym = enzymUsed, missedCleavages = missedCleavages, unique = TRUE) %>% 
    # Unlist to generate the list we can use to look-up fragments.
    BiocGenerics::unlist() %>% BiocGenerics::unique() %>% Biostrings::sort()
  
  
  # Compare against TxDb and additional peptides. ---------------------------
  
  ParallelLogger::logTrace('\tProteoDisco - Checking unique fragments (against reference database(s))')
  
  # Helper function to check number of unique fragments and which fragments those are.
  predictNumberOfUniqueFragments <- function(data, enzymUsed, missedCleavages, TxDb.Fragments, checkWithinMutantSeqs){
    
    # Remove previously-annotated columns.
    data$proteotypicFragments <- data$proteotypicFragmentsCount <- data$totalFragments <- NULL
    
    # Add a look-up row-index.
    data$rowIndex <- base::seq_len(base::nrow(data))
    seqMut <- data$AA.SequenceMut
    base::names(seqMut) <- data$rowIndex
    
    mut.Fragments <- cleaver::cleave(seqMut, enzym = enzymUsed, missedCleavages = missedCleavages, unique = TRUE) %>% 
      BiocGenerics::unlist() %>% 
      Biostrings::sort()
    
    # Count the number of (cleaved) fragments per row.
    resultCleaving <- tibble::as_tibble(base::as.data.frame(base::table(base::names(mut.Fragments)), responseName = 'totalFragments'))
    
    # Unique mutant-fragments should only match with themselves (if they are unique).
    if(checkWithinMutantSeqs){
      mut.Fragments <- mut.Fragments[S4Vectors::countMatches(mut.Fragments, mut.Fragments) == 1,]
    }
    
    # Retrieve the proteotypic fragments not present in the reference database.
    mut.Fragments.Proteotypic <- base::unique(mut.Fragments[base::is.na(BiocGenerics::match(mut.Fragments, TxDb.Fragments, nomatch = NA_integer_))])
    S4Vectors::mcols(mut.Fragments.Proteotypic)$rowIndex <- base::names(mut.Fragments.Proteotypic)
    
    # Melt to rowIndex:proteotypicFragment(s) and count nr. of proteotypic fragments.
    mut.Fragments.Proteotypic <- tibble::as_tibble(S4Vectors::DataFrame(mut.Fragments.Proteotypic)) %>% 
      dplyr::group_by(.data$rowIndex) %>% 
      dplyr::summarise(
        proteotypicFragments = paste(base::unique(.data$X), collapse = ', '),
        proteotypicFragmentsCount = dplyr::n_distinct(.data$X)
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::full_join(resultCleaving, by = c('rowIndex' = 'Var1'))
    
    # Merge results back to original data.
    dataWithFragments <- S4Vectors::merge(data, S4Vectors::DataFrame(mut.Fragments.Proteotypic), by = 'rowIndex', all.x = TRUE)
    dataWithFragments$Tx.SequenceMut <- Biostrings::DNAStringSet(dataWithFragments$Tx.SequenceMut)
    dataWithFragments$AA.SequenceMut <- Biostrings::AAStringSet(dataWithFragments$AA.SequenceMut)
    
    # Return.
    return(dataWithFragments)
    
  }
  
  # Gather all mutant transcripts.
  transcripts <- ProteoDisco::mutantTranscripts(x)
  
  # Check imported genomic variants.
  if(base::nrow(transcripts$genomicVariants) > 0){
    transcripts$genomicVariants <- predictNumberOfUniqueFragments(transcripts$genomicVariants, enzymUsed, missedCleavages, TxDb.Fragments, checkWithinMutantSeqs)
    x <- setMutantTranscripts(x, transcripts = transcripts$genomicVariants, slotType = 'genomicVariants')  
  }
  
  # Check manually-inputted sequences.
  if(base::nrow(transcripts$manualSequences) > 0){
    transcripts$manualSequences <- predictNumberOfUniqueFragments(transcripts$manualSequences, enzymUsed, missedCleavages, TxDb.Fragments, checkWithinMutantSeqs)
    x <- setMutantTranscripts(x, transcripts = transcripts$manualSequences, slotType = 'manualSequences')
  }
  
  
  # Return ------------------------------------------------------------------
  
  return(x)
})
