context("Determination and incorporation of proteotypic fragments")

test_that("Proteotypic fragments can be determined and incorporated", {
  
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
    ignoreNonMatch = TRUE
  )
  
  # Incorporate genomic variants.
  ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
    ProteoDiscography = ProteoDiscography.hg19,
    aggregateSamples = FALSE, 
    aggregateWithinExon = TRUE, 
    aggregateWithinTranscript = TRUE, 
    ignoreOverlappingMutations = TRUE, 
    threads = 1
  )
  
  # Check proteotypic fragments.
  ProteoDiscography.hg19 <- ProteoDisco::checkProteotypicFragments(
    ProteoDiscography.hg19, 
    checkWithinMutantSeqs = TRUE
  )
  
  # Check correct nr. of proteotypic fragments.
  testthat::expect_equal(
    object = sum(ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$genomicVariants$proteotypicFragmentsCount > 0, na.rm = TRUE),
    expected = 20
  )
  
})
