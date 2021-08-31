context("Importing manual sequences as additional transcripts")

test_that("Manual transcripts are incorporated as expected", {

  # Import example ProteoDiscography (hg19)
  data('ProteoDiscographyExample.hg19', package = 'ProteoDisco')
  ProteoDiscographyExample.hg19 <- ProteoDisco::setTxDb(ProteoDiscographyExample.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  ProteoDiscographyExample.hg19 <- ProteoDisco::setGenomicSequences(ProteoDiscographyExample.hg19, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

  # Retrieve TMPRSS2-ERG sequence from ENA.
  testSeq1 <- Biostrings::DNAString('ATGACCGCGTCCTCCTCCAGCGACTATGGACAGACTTCCAAGATGAGCCCACGCGTCCCTCAGCAGGATTGGCTGTCTCAACCCCCAGCCAGGGTCACCATCAAAATGGAATGTAACCCTAGCCAGGTGAATGGCTCAAGGAACTCTCCTGATGAATGCAGTGTGGCCAAAGGCGGGAAGATGGTGGGCAGCCCAGACACCGTTGGGATGAACTACGGCAGCTACATGGAGGAGAAGCACATGCCACCCCCAAACATGACCACGAACGAGCGCAGAGTTATCGTGCCAGCAGATCCTACGCTATGGAGTACAGACCATGTGCGGCAGTGGCTGGAGTGGGCGGTGAAAGAATATGGCCTTCCAGACGTCAACATCTTGTTATTCCAGAACATCGATGGGAAGGAACTGTGCAAGATGACCAAGGACGACTTCCAGAGGCTCACCCCCAGCTACAACGCCGACATCCTTCTCTCACATCTCCACTACCTCAGAGAGACTCCTCTTCCACATTTGACTTCAGATGATGTTGATAAAGCCTTACAAAACTCTCCACGGTTAATGCATGCTAGAAACACAGGGGGTGCAGCTTTTATTTTCCCAAATACTTCAGTATATCCTGAAGCTACGCAAAGAATTACAACTAGGCCAGTCTCTTACAGATAA')

  # Retrieve partial CDS of BCR-ABL from ENA.
  testSeq2 <- Biostrings::DNAString('ATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAGGTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTTTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAAT')

  # Generate DataFrame containing the mutant sequences and metadata.
  manualSeq <- S4Vectors::DataFrame(
    sample = rep('example', 2),
    identifier = c('TMPRSS2-ERG prostate cancer specific isoform 1', 'partial bcr/abl e8a2 fusion CDS'),
    gene = c('TMPRSS2-ERG', 'BCR-ABL'),
    Tx.SequenceMut = Biostrings::DNAStringSet(base::list(testSeq1, testSeq2))
  )

  ProteoDiscographyExample.hg19 <- ProteoDisco::importTranscriptSequences(
    ProteoDiscographyExample.hg19,
    transcripts = manualSeq,
    removeExisting = TRUE
  )

  testthat::expect_equal(
    object = nrow(ProteoDisco::getDiscography(ProteoDiscographyExample.hg19)$manualSequences),
    expected = 2
  )

  testthat::expect_equal(
    object = nrow(ProteoDisco::mutantTranscripts(ProteoDiscographyExample.hg19)$manualSequences),
    expected = 2
  )

})
