#' @title Display a summary of ProteoDiscography
#' @param object (ProteoDiscography): ProteoDiscography object.
#'
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' print(ProteoDiscography.hg19)
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
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' show(ProteoDiscography.hg19)
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
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' summary(ProteoDiscography.hg19)
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
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' seqinfo(ProteoDiscography.hg19)
#' 
#' @return Returns an Seqinfo-class containing information on the imported genomic sequences.
#' @rdname seqinfo
#' @export
setGeneric('seqinfo', function(x) standardGeneric("seqinfo"))

#' @title Retrieve seqlevels.
#' @param x (ProteoDiscography): ProteoDiscography object.
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' seqlevels(ProteoDiscography.hg19)
#' 
#' @return seqlevels of the improted genomic sequences.
#' @rdname seqlevels
#' @export
setGeneric('seqlevels', function(x) standardGeneric("seqlevels"))

#' @title Display the organism of the TxDb
#' @param object (ProteoDiscography): ProteoDiscography object.
#' @examples
#' 
#' ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' )
#' 
#' organism(ProteoDiscography.hg19)
#' 
#' @return Returns the specified organism during creation.
#' @rdname organism
#' @export
setGeneric('organism', function(object) standardGeneric("organism"))

#' @title Retrieve imported genomic variants, splice-junctions and manual sequences.
#' @param x (ProteoDiscography): ProteoDiscography object.
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
#' # Retrieve the imported records.
#' getDiscography(ProteoDiscography.hg19)
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
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography.hg19)$genomicVariants[1:10], slotType = 'genomicVariants')
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography.hg19)$spliceJunctions[1:10], slotType = 'spliceJunctions')
#' # ProteoDisco::setMutantTranscripts(ProteoDiscography.hg19)$manualSequences[1:10], slotType = 'manualSequences')
#' 
#' # Example using genomic variants
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
#' # Only keep the first10 records.
#' ProteoDiscography.hg19 <- ProteoDisco::setMutantTranscripts(
#'   x = ProteoDiscography.hg19, 
#'   transcripts = ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$genomicVariants[1:10,], 
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
#' # Retrieve all generated mutant transcripts.
#' mutantTranscripts(ProteoDiscography.hg19)
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
#'  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
#'    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'  )
#'  
#'  ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
#'    ProteoDiscography = ProteoDiscography.hg19,
#'    # Provide the VCF / MAF files, if more then one supply a vector of files and corresponding samplenames.
#'    files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'), 
#'    # We can replace the original samples within the VCF with nicer names.
#'    samplenames = 'Validation Set (GRCh37)',
#'    # Number of threads used for parallelization.
#'    # We run samples sequentially and parallelize within (variant-wise multi-threading).
#'    threads = 1, 
#'    # To increase import-speed for this example, do not check for validity of the reference anchor with the given reference sequences.
#'    performAnchorCheck = FALSE
#'  )
#'  
#'  # Results will now contain additional information about proteotypic fragments.
#'  ProteoDiscography.hg19 <- checkProteotypicFragments(ProteoDiscography.hg19)
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
