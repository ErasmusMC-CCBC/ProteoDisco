context("Importing genomic SNV, InDel and MNV from VCF")

test_that("Variants in VCF file are incorporated as expected", {
  
  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  )
  
  # One or more VCF or MAF files.
  files.muts <- system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco')
  
  # Import all mutations from one or multiple VCF (or MAF) files.
  ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
    ProteoDiscography.hg19, 
    files.muts, 
    ignoreNonMatch = TRUE,
    threads = 1
  )
  
  ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
    ProteoDiscography = ProteoDiscography.hg19,
    aggregateSamples = FALSE, 
    aggregateWithinExon = TRUE, 
    aggregateWithinTranscript = TRUE, 
    ignoreOverlappingMutations = TRUE, 
    threads = 1
  )
  
  testSetProteins_hg19 <- readr::read_csv(file = system.file("extdata", "Table1_testSetProteins.csv", package = "ProteoDisco")) %>% dplyr::filter(Build == 'hg19')
  
  # Missing two out of 28 at this point, as expected
  testthat::expect_equal(sum(!is.na(BiocGenerics::match(testSetProteins_hg19$`Mutant Protein`, ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$genomicVariants$AA.SequenceMut))), expected = 26)
  
  # Add the missing two sequences.
  ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
    ProteoDiscography = ProteoDiscography.hg19,
    aggregateSamples = FALSE,
    aggregateWithinExon = FALSE,
    aggregateWithinTranscript = TRUE,
    ignoreOverlappingMutations = TRUE,
    threads = 1
  )
  
  # All 28 should be found.
  testthat::expect_equal(sum(!is.na(BiocGenerics::match(testSetProteins_hg19$`Mutant Protein`, ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$genomicVariants$AA.SequenceMut))), expected = 28)
  
})

test_that("SJs are incorporated as expected", {
  
  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  )
  
  files.Splicing <- c(
    system.file('extdata', 'spliceJunctions_pyQUILTS_chr22.bed', package = 'ProteoDisco')
  )
  
  # Import splice-junctions from BED or SJ,out.tab files into our ProteoDiscography.
  ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
    ProteoDiscography = ProteoDiscography.hg19, 
    inputSpliceJunctions = files.Splicing, 
    # (Optional) Rename samples.
    samples = 'pyQUILTS',
    isTopHat = TRUE,
    aggregateSamples = FALSE,
    removeExisting = TRUE
  )
  
  # Or, import splice-junctions (even spanning different chromosomes) based on our format.
  testSJ <- readr::read_tsv(system.file('extdata', 'validationSetSJ_hg19.txt', package = 'ProteoDisco')) %>% 
    dplyr::select(sample, identifier, junctionA, junctionB)
  
  # Add custom SJ to ProteoDiscography.
  ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
    ProteoDiscography = ProteoDiscography.hg19, 
    inputSpliceJunctions = testSJ, 
    # Append to existing SJ-input.
    removeExisting = FALSE
  )
  
  # Generate junction-models from non-canonical splice-junctions.
  ProteoDiscography.hg19 <- ProteoDisco::generateJunctionModels(
    ProteoDiscography = ProteoDiscography.hg19, 
    maxDistance = 150, 
    skipCanonical = TRUE,
    threads = 1,
  )
  
  testthat::expect_equal(length(ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$spliceJunctions$identifier), 70, info = "Should have 70 non-canonical splice-transcripts from the two sources")
  
})
