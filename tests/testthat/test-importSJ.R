context("Importing SJ")

test_that("SJs are incorporated as expected", {

  # Import example ProteoDiscography (hg19) and re-link TxDb.
  data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
  ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

  files.Splicing <- c(
    system.file('extdata', 'spliceJunctions_pyQUILTS_chr22.bed', package = 'ProteoDisco')
  )

  # Import splice-junctions from BED or SJ,out.tab files into our ProteoDiscography.
  ProteoDiscographyExample.hg19 <- ProteoDisco::importSpliceJunctions(
    ProteoDiscography = ProteoDiscographyExample.hg19,
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
  ProteoDiscographyExample.hg19 <- ProteoDisco::importSpliceJunctions(
    ProteoDiscography = ProteoDiscographyExample.hg19,
    inputSpliceJunctions = testSJ,
    # Append to existing SJ-input.
    removeExisting = FALSE
  )

  # Generate junction-models from non-canonical splice-junctions.
  ProteoDiscographyExample.hg19 <- ProteoDisco::generateJunctionModels(
    ProteoDiscography = ProteoDiscographyExample.hg19,
    maxDistance = 150,
    skipCanonical = TRUE,
    threads = 1
  )

  testthat::expect_equal(length(ProteoDisco::mutantTranscripts(ProteoDiscographyExample.hg19)$spliceJunctions$identifier), 70, info = "Should have 70 non-canonical splice-transcripts from the two sources")

})
