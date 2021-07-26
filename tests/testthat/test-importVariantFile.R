context("Importing genomic SNV, InDel and MNV from VCF")

test_that("VCF file is placed in correct ProteoDiscography slot", {
  
  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  )
  
  # One or more VCF or MAF files.
  files.muts <- system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco')
  
  # Import all mutations from one or multiple VCF (or MAF) files.
  ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
    ProteoDiscography.hg19, 
    files.muts, ignoreNonMatch = FALSE
  )
  
  z <- ProteoDisco::summary(ProteoDiscography.hg19, verbose = FALSE)
  
  testthat::expect_is(ProteoDiscography.hg19, "ProteoDiscography", info = "Should be ProteoDiscography object.")
  testthat::expect_equal(z$overviewMutations$InDel, 10, info = "Should have 10 accessible InDel variants.")
  testthat::expect_equal(z$overviewMutations$MNV, 1, info = "Should have 1 accessible MNV variants.")
  testthat::expect_equal(z$overviewMutations$SNV, 17, info = "Should have 17 accessible SNV variants.")
  
})

test_that("Duplicate sample is not allowed (or overwritten)", {
  
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
    overwriteDuplicateSamples = FALSE
  )
  
  #Expect an error to pop due to duplicate samples and overwriteDuplicateSamples set to FALSE
  testthat::expect_error(ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
    ProteoDiscography.hg19, 
    files.muts, 
    ignoreNonMatch = TRUE, 
    overwriteDuplicateSamples = FALSE)
  )
  
  #Expect no error pop due to duplicate samples and overwriteDuplicateSamples set to TRUE
  testthat::expect_error(ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
    ProteoDiscography.hg19, 
    files.muts,
    ignoreNonMatch = TRUE, 
    overwriteDuplicateSamples = TRUE),
    NA
  )
  
})
