context("Converting TxID to gene SYMBOL")

test_that("TxIDs are correctly converted to gene SYMBOLs", {
  
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
  
  # Retrieve all imported mutant transcripts.
  transcripts <- ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)
  testthat::expect_null(object = transcripts$genomicVariants$SYMBOL)
  
  # Convert ENTREZ identifiers to gene symbols.
  geneSymbols <- unique(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = transcripts$genomicVariants$geneID, 
    keytype = 'ENTREZID', columns = 'SYMBOL')
  )
  
  # Add the gene symbols back to the mutant transcript.
  transcripts$genomicVariants <- S4Vectors::DataFrame(base::merge(transcripts$genomicVariants, geneSymbols, by.x = 'geneID', by.y = 'ENTREZID', all.x = TRUE))
  
  # Add the mutant transcripts (with symbols) back into the ProteoDiscography.
  ProteoDiscography.hg19 <- ProteoDisco::setMutantTranscripts(
    ProteoDiscography.hg19, 
    transcripts = transcripts$genomicVariants, 
    slotType = 'genomicVariants'
  )
  
  testthat::expect_equal(
    object = length(unique(transcripts$genomicVariants$SYMBOL)), 
    expected = 17
  )
  
})
