context("Importing genomic SNV, InDel and MNV from VCF")

test_that("Min. import to export workflow", {


  # Generate ProteoDiscography ----------------------------------------------

  ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  )


  # Import validation set ---------------------------------------------------

  files.muts <- system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco')

  # Import all mutations from one or multiple VCF (or MAF) files.
  ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
    ProteoDiscography.hg19,
    files.muts,
    ignoreNonMatch = FALSE,
    threads = 1
  )

  z <- ProteoDisco::summary(ProteoDiscography.hg19, verbose = FALSE)

  testthat::expect_is(ProteoDiscography.hg19, "ProteoDiscography", info = "Should be ProteoDiscography object.")
  testthat::expect_equal(z$overviewMutations$InDel, 10, info = "Should have 10 accessible InDel variants.")
  testthat::expect_equal(z$overviewMutations$MNV, 1, info = "Should have 1 accessible MNV variants.")
  testthat::expect_equal(z$overviewMutations$SNV, 17, info = "Should have 17 accessible SNV variants.")


  # Incorporate genomic variants --------------------------------------------

  testSetProteins_hg19 <- readr::read_csv(file = system.file("extdata", "Table1_testSetProteins.csv", package = "ProteoDisco")) %>% dplyr::filter(Build == 'hg19')

  ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
    ProteoDiscography = ProteoDiscography.hg19,
    aggregateSamples = FALSE,
    aggregateWithinExon = TRUE,
    aggregateWithinTranscript = TRUE,
    ignoreOverlappingMutations = TRUE,
    threads = 1
  )

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

  # Try export to FASTA.
  ProteoDisco::exportProteoDiscography(ProteoDiscography.hg19, outFile = "test.FASTA")
  testthat::expect_equal(object = list.files(pattern = "test.FASTA"), expected = "test.FASTA")
  file.remove("test.FASTA")
  testthat::expect_equal(object = list.files(pattern = "test.FASTA"), expected = character(0))

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
