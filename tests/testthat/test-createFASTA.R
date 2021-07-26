context("Create FASTA file")

test_that("FASTA file can be created", {
  
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
  
  # Check proteotypic fragments.
  ProteoDiscography.hg19 <- ProteoDisco::checkProteotypicFragments(
    ProteoDiscography.hg19, 
    checkWithinMutantSeqs = FALSE
  )
  
  testthat::expect_equal(
    object = length(which(ProteoDiscography.hg19@mutantTranscripts.genomicVariants$proteotypicFragmentsCount > 0)), 
    expected = 32
  )
  
  ProteoDiscography.hg19 <- ProteoDisco::checkProteotypicFragments(
    ProteoDiscography.hg19, 
    checkWithinMutantSeqs = TRUE
  )
  
  testthat::expect_equal(
    object = length(which(ProteoDiscography.hg19@mutantTranscripts.genomicVariants$proteotypicFragmentsCount > 0)), 
    expected = 20
  )
  
  ProteoDisco::exportProteoDiscography(ProteoDiscography.hg19, outFile = "test.FASTA", minProteotypicFragments = 1)
  testthat::expect_equal(object = list.files(pattern = "test.FASTA"), expected = "test.FASTA")
  file.remove("test.FASTA")
  testthat::expect_equal(object = list.files(pattern = "test.FASTA"), expected = character(0))
})
