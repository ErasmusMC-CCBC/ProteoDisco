context("Determination and incorporation of proteotypic fragments")

test_that("Proteotypic fragments can be determined and incorporated", {

    # Import example ProteoDiscography (hg19) and re-link TxDb.
    data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
    ProteoDiscographyExample.hg19@TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

    # Check proteotypic fragments.
    ProteoDiscographyExample.hg19 <- ProteoDisco::checkProteotypicFragments(
        ProteoDiscographyExample.hg19,
        checkWithinMutantSeqs = TRUE
    )

    # Check correct nr. of proteotypic fragments.
    testthat::expect_equal(
        object = sum(ProteoDisco::mutantTranscripts(ProteoDiscographyExample.hg19)$genomicVariants$proteotypicFragmentsCount > 0, na.rm = TRUE),
        expected = 20
    )

})
