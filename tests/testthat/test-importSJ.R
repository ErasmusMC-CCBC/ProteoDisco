context("Importing SJ")

test_that("SJs are incorporated as expected", {

  # Import example ProteoDiscography (hg19)
  data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
  ProteoDiscographyExample.hg19 <- ProteoDisco::setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  ProteoDiscographyExample.hg19 <- ProteoDisco::setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

  # Import splice-junctions (even spanning different chromosomes) based on our format.
  testSJ <- readr::read_tsv(system.file('extdata', 'validationSetSJ_hg19.txt', package = 'ProteoDisco')) %>%
    dplyr::select(sample, identifier, junctionA, junctionB)

  # Add custom SJ to ProteoDiscography.
  ProteoDiscographyExample.hg19 <- ProteoDisco::importSpliceJunctions(
    ProteoDiscography = ProteoDiscographyExample.hg19,
    inputSpliceJunctions = testSJ,
    # Append to existing SJ-input.
    removeExisting = TRUE
  )

  # Generate junction-models from non-canonical splice-junctions.
  ProteoDiscographyExample.hg19 <- ProteoDisco::generateJunctionModels(
    ProteoDiscography = ProteoDiscographyExample.hg19,
    maxDistance = 150,
    skipCanonical = TRUE,
    threads = 1
  )

  testthat::expect_equal(length(ProteoDisco::mutantTranscripts(ProteoDiscographyExample.hg19)$spliceJunctions$identifier), 5, info = "Should have 5 non-canonical splice-transcripts.")

})
