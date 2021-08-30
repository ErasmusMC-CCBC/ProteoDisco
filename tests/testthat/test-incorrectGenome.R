context("Importing variants from VCF-file with incorrect genome build.")

test_that("Non-matching reference anchors", {

  ProteoDiscography.hg38 <- ProteoDisco::generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )

  # Attempt import, but break on mismatches with reference genome
  testthat::expect_error(
    ProteoDisco::importGenomicVariants(
      ProteoDiscography.hg38,
      system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'),
      ignoreNonMatch = FALSE
    )
  )

  # Unless ignore is set to TRUE.
  testthat::expect_s4_class(
    ProteoDisco::importGenomicVariants(
      ProteoDiscography.hg38,
      system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'),
      ignoreNonMatch = TRUE
    ), class = 'ProteoDiscography'
  )

})
