## ----knitr_init, echo = FALSE, cache = FALSE, results = 'hide', warning=FALSE, message=FALSE----------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(ProteoDisco))

## ---- label = 'Minimal Workflow', eval = TRUE, results = 'hide', message = FALSE----------------------------------------------------------------------------------------------------

# Generate the ProteoDiscography (containing the genome-sequences and annotations used throughout the incorporation of genomic variants)
# In this case, we supply an existing reference genome with known annotations (GRCh37 with GENCODE annotations).
ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
  genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
)

# Supply the ProteoDiscography with genomic variants to incorporate in downstream analysis. This can be one or multiple VCF / MAF files.
# Additional manual sequences and exon-exon mapping (i.e., splice junctions) can also be given as shown in the sections below.
ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
  ProteoDiscography = ProteoDiscography.hg19,
  # Provide the VCF / MAF files.
  files = system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco'), 
  # We can replace the original samples within the VCF with nicer names.
  samplenames = 'Validation Set (GRCh37)',
  # Number of threads used for parallelization.
  # We run samples sequentially and parallelize within (variant-wise multi-threading).
  threads = 1, 
  # To increase import-speed for this example, do not check for validity of the reference anchor with the given reference sequences.
  performAnchorCheck = FALSE
)

# Incorporate the given genomic variants into their respective overlapping coding-sequences (i.e. transcripts).
# This can be done in a sample-specific or aggregated cohort-wide manner and can be performed per exon or transcript separately.
ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
  ProteoDiscography = ProteoDiscography.hg19,
  # Do not aggregate samples and generate mutant transcripts from the mutations per sample.
  aggregateSamples = FALSE, 
  # If there are multiple mutations within the same exon (CDS), place them on the same mutant CDS sequence.
  aggregateWithinExon = TRUE, 
  # Aggregate multiple mutant exons (CDS) within the same transcripts instead of incorporating one at a time.
  aggregateWithinTranscript = TRUE, 
  # If there are overlapping mutations on the same coding position, retain only the first of the overlapping mutations.
  # If set to FALSE, throw an error and specify which CDS had overlapping mutations.
  ignoreOverlappingMutations = TRUE, 
  # Number of threads.
  threads = 1
)

## ---- label = 'Export FASTA', eval = FALSE, results = 'hide'------------------------------------------------------------------------------------------------------------------------
#  # Output custom protein database (FASTA) containing the generated protein variant sequences.
#  ProteoDisco::exportProteoDiscography(
#    ProteoDiscography = ProteoDiscography.hg19,
#    outFile = 'example.fasta',
#    aggregateSamples = TRUE
#  )

## ---- label = 'Logging options', eval = TRUE, results = 'hide'----------------------------------------------------------------------------------------------------------------------
# Display tracing messages:
ParallelLogger::clearLoggers()
ParallelLogger::registerLogger(ParallelLogger::createLogger(threshold = 'TRACE', appenders =list(ParallelLogger::createConsoleAppender(layout = ParallelLogger::layoutTimestamp))))

# Or only display general information messages:
ParallelLogger::clearLoggers()
ParallelLogger::registerLogger(ParallelLogger::createLogger(threshold = 'INFO', appenders =list(ParallelLogger::createConsoleAppender(layout = ParallelLogger::layoutTimestamp))))

## ---- label = 'TxDb from FASTA and GTF', eval = FALSE-------------------------------------------------------------------------------------------------------------------------------
#  
#  # Download the latest annotation.
#  URL.annotations <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gff3.gz'
#  utils::download.file(URL.annotations, basename(URL.annotations))
#  
#  # Download the respective reference genome (GRCh38.p12)
#  URL.genome <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p12.genome.fa.gz'
#  utils::download.file(URL.genome, basename(URL.genome))
#  
#  # Use the latest annotations to generate a TxDb.
#  TxDb <- GenomicFeatures::makeTxDbFromGFF(
#    file = basename(URL.annotations),
#    dataSource = 'GENCODE (v37)',
#    organism = 'Homo sapiens', taxonomyId = 9606
#  )
#  
#  # Import the genome-sequences.
#  genomeSeqs <- Biostrings::readDNAStringSet(filepath = basename(URL.genome), format = 'fasta')
#  
#  # Clean the chromosomal names and use the chr-prefix.
#  names(genomeSeqs) <- gsub(' .*', '', names(genomeSeqs))
#  
#  # Generate the ProteoDiscography.
#  ProteoDiscography <- ProteoDisco::generateProteoDiscography(
#    TxDb = TxDb,
#    genomeSeqs = genomeSeqs,
#    geneticCode = 'Standard'
#  )

## ---- label = 'Generate ProteoDiscrography from GRCh37 resources', eval = TRUE, results = 'hide'------------------------------------------------------------------------------------
# Attach the required packages.
require(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Generate the ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
  genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
)

## ---- label = 'Inspection', eval = TRUE, results = 'hide'---------------------------------------------------------------------------------------------------------------------------
ProteoDiscography.hg19 # or print(ProteoDiscography.hg19)
ProteoDisco::summary(ProteoDiscography.hg19)

## ---- label = 'Generate ProteoDiscrography (GRCh37) for subsequent testing', eval = TRUE, results = 'hide'--------------------------------------------------------------------------
# Lets produce a ProteoDiscography using pre-generated (BioConductor) resources for hg19 with UCSC annotations.
require(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Initialize the ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::generateProteoDiscography(
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
  genomeSeqs = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
)

# Select one or more VCF or MAF files.
# In this case, we make use of the provided GRCh37 test file.
files.muts <- c(
  system.file('extdata', 'validationSet_hg19.vcf', package = 'ProteoDisco')
)

# Import the genomic (somatic) variants.
ProteoDiscography.hg19 <- ProteoDisco::importGenomicVariants(
  ProteoDiscography = ProteoDiscography.hg19,
  # Provide the VCF / MAF files.
  files = files.muts, 
  # We can replace the original samples within the VCF with nicer names.
  samplenames = 'Validation Set (GRCh37)',
  # For illustration purposes, do not check the validity of the reference anchor.
  performAnchorCheck = FALSE
)


## ---- label = 'Inspect example ProteoDiscrography (GRCh37)', eval = TRUE, results = 'hide'------------------------------------------------------------------------------------------
library(ggplot2)

# Print an overview of supplied variants, their respective mutational types and an overview per patients.
ProteoDisco::summary(ProteoDiscography.hg19, verbose = FALSE)

# We can also capture the summary output and display / handle only certain parts. 
# The 'verbose' parameters toggles whether additional information is printed to the console.
z <- ProteoDisco::summary(ProteoDiscography.hg19, verbose = FALSE)

# Generate a quick overview of the mutational categories (SNV, MNV and InDels) per imported sample.
reshape2::melt(z$overviewMutations, id.vars = 'sample') %>% 
  dplyr::filter(!is.na(value)) %>% 
  ggplot2::ggplot(ggplot2::aes(x = sample, y = value, fill = variable)) + 
  ggplot2::geom_bar(stat = 'identity', width = .8, color = 'black', position = ggplot2::position_dodge(preserve = 'single')) + 
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 10000)) +
  ggplot2::labs(x = 'Provided samples', y = 'Nr. of variants (log10)') +
  ggplot2::scale_fill_brewer(name = 'Category', palette = 'Greens') +
  ggplot2::theme_minimal()

## ---- label = 'Add manual transcript sequences', eval = TRUE, results = 'hide', cache = TRUE----------------------------------------------------------------------------------------

# TMPRSS2-ERG sequence from ENA.
testSeq1 <- Biostrings::DNAString('ATGACCGCGTCCTCCTCCAGCGACTATGGACAGACTTCCAAGATGAGCCCACGCGTCCCTCAGCAGGATTGGCTGTCTCAACCCCCAGCCAGGGTCACCATCAAAATGGAATGTAACCCTAGCCAGGTGAATGGCTCAAGGAACTCTCCTGATGAATGCAGTGTGGCCAAAGGCGGGAAGATGGTGGGCAGCCCAGACACCGTTGGGATGAACTACGGCAGCTACATGGAGGAGAAGCACATGCCACCCCCAAACATGACCACGAACGAGCGCAGAGTTATCGTGCCAGCAGATCCTACGCTATGGAGTACAGACCATGTGCGGCAGTGGCTGGAGTGGGCGGTGAAAGAATATGGCCTTCCAGACGTCAACATCTTGTTATTCCAGAACATCGATGGGAAGGAACTGTGCAAGATGACCAAGGACGACTTCCAGAGGCTCACCCCCAGCTACAACGCCGACATCCTTCTCTCACATCTCCACTACCTCAGAGAGACTCCTCTTCCACATTTGACTTCAGATGATGTTGATAAAGCCTTACAAAACTCTCCACGGTTAATGCATGCTAGAAACACAGGGGGTGCAGCTTTTATTTTCCCAAATACTTCAGTATATCCTGAAGCTACGCAAAGAATTACAACTAGGCCAGTCTCTTACAGATAA')

# Partial CDS of BCR-ABL from ENA.
testSeq2 <- Biostrings::DNAString('ATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAGGTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTTTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAAT')

# Generate DataFrame containing the mutant sequences and metadata and add these to our test sample.
manualSeq <- S4Vectors::DataFrame(
  sample = rep('Validation Set (GRCh37)', 2),
  identifier = c('TMPRSS2-ERG prostate cancer specific isoform 1', 'Homo sapiens (human) partial bcr/c-abl oncogene protein'),
  gene = c('TMPRSS2-ERG', 'BCR-ABL1'),
  Tx.SequenceMut = Biostrings::DNAStringSet(base::list(testSeq1, testSeq2))
)

# Add to ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::importTranscriptSequences(
  ProteoDiscography.hg19, 
  transcripts = manualSeq,
  removeExisting = TRUE
)

## ---- label = 'View imported ProteoDiscography records', eval = TRUE, results = 'hide'----------------------------------------------------------------------------------------------

# Retrieve / view the records imported into the ProteoDiscography.
getDiscography(ProteoDiscography.hg19)

## ---- label = 'Incorporate genomic variants into ProteoDiscrography (GRCh37)', eval = TRUE, results = 'hide', message = FALSE-------------------------------------------------------

ProteoDiscography.hg19 <- ProteoDisco::incorporateGenomicVariants(
  ProteoDiscography = ProteoDiscography.hg19,
  # Do not aggregate samples and generate mutant transcripts from the mutations per sample.
  aggregateSamples = FALSE, 
  # If there are multiple mutations within the same exon (CDS), place them on the same mutant CDS sequence.
  aggregateWithinExon = TRUE, 
  # Aggregate multiple mutant exons (CDS) within the same transcripts instead of incorporating one at a time.
  aggregateWithinTranscript = TRUE, 
  # If there are overlapping mutations on the same coding position, retain only the first of the overlapping mutations.
  # If set to FALSE, throw an error and specify which CDS had overlapping mutations.
  ignoreOverlappingMutations = TRUE, 
  # Number of threads.
  threads = 4
)


## ---- label = 'Inspection post-incorporation', eval = TRUE, results = 'hide'--------------------------------------------------------------------------------------------------------

# Simply print the ProteoDiscography and it will state the number of mutant transcripts currently generated.
ProteoDisco::summary(ProteoDiscography.hg19)

# Retrieve the mutant transcript sequences.
# This will retrieve a list of all the generated (and imported) transcript sequences per input type and the respective amino-acid sequence.
x <- ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)

# Visualization of the number of incorporated genomic mutations per transcript.
barplot(table(x$genomicVariants$numberOfMutationsInTx))


## ---- label = 'Subsetting', eval = TRUE, results = 'hide'---------------------------------------------------------------------------------------------------------------------------
x <- ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)
x <- x$genomicVariants
x

# Only keep the first 10 mutant transcripts of the incorporated genomic variants.
ProteoDiscography.hg19.subset <- ProteoDisco::setMutantTranscripts(ProteoDiscography.hg19, x[1:10,], slotType = 'genomicVariants')

y <- ProteoDisco::mutantTranscripts(ProteoDiscography.hg19.subset)
y <- y$genomicVariants
y


## ---- label = 'Import splice-junctions', eval = TRUE, results = 'hide', message = FALSE---------------------------------------------------------------------------------------------
# Import splice-junctions from BED files.
files.Splicing <- c(
  system.file('extdata', 'spliceJunctions_pyQUILTS_chr22.bed', package = 'ProteoDisco')
)

# Import splice-junctions from BED or SJ,out.tab files into our ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
  ProteoDiscography = ProteoDiscography.hg19, 
  inputSpliceJunctions = system.file('extdata', 'spliceJunctions_pyQUILTS_chr22.bed', package = 'ProteoDisco'), 
  # (Optional) Rename samples.
  samples = c('pyQUILTS'), 
  # Specify that the given BED files are obtained from TopHat.
  # Chromosomal coordinates from TopHat require additional formatting.
  isTopHat = TRUE,
  aggregateSamples = FALSE,
  removeExisting = TRUE
)

# Or, import splice-junctions (even spanning different chromosomes) based on our format.
testSJ <- readr::read_tsv(system.file('extdata', 'validationSetSJ_hg19.txt', package = 'ProteoDisco')) %>% 
  dplyr::select(sample, identifier, junctionA, junctionB)

# Add custom SJ to ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::importSpliceJunctions(
  ProteoDiscography = ProteoDiscography.hg19, 
  inputSpliceJunctions = testSJ, 
  # Remove existing SJ-input.
  removeExisting = TRUE
)

# Generate junction-models from non-canonical splice-junctions.
ProteoDiscography.hg19 <- ProteoDisco::generateJunctionModels(
  ProteoDiscography = ProteoDiscography.hg19, 
  # Max. distance from a known exon-boundary before introducing a novel exon.
  # If an adjacent exon is found within this distance, it will shorten or elongate that exon towards the SJ.
  maxDistance = 150, 
  # Should we skip known exon-exon junctions (in which both the acceptor and donor are located on known adjacent exons within the same transcript)
  skipCanonical = TRUE,
  # Perform on multiple threads (optional)
  threads = 4
)

## ---- label = 'Check proteotypic fragments.', eval = FALSE, results = 'hide'--------------------------------------------------------------------------------------------------------
#  # Load additional peptide sequences against which proteotypic fragments are checked.
#  # In this case, load the first 100 human Swiss-Prot/Uniprot database entries.
#  peptides.SwissProt <- base::textConnection(RCurl::getURL('https://www.uniprot.org/uniprot/?query=*&format=fasta&limit=100&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes')) %>%
#    seqinr::read.fasta(., seqtype = 'AA', seqonly = FALSE, as.string = TRUE) %>%
#    Biostrings::AAStringSetList() %>%
#    unlist()
#  
#  # Check proteotypic fragments.
#  ProteoDiscography.hg19 <- ProteoDisco::checkProteotypicFragments(
#    ProteoDiscography.hg19,
#    additionalPeptides = peptides.SwissProt,
#    # Protease used in the experiment (see cleaver package)
#    enzymUsed = 'trypsin',
#    # Number of potentially-missed cleavages (see cleaver package)
#    missedCleavages = 1,
#    # Should proteotypic fragments also be unique among the generated variant-sequences.
#    checkWithinMutantSeqs = FALSE
#  )
#  
#  # We can now appreciate that additional information is added to the mutant transcript denoting how many, and which fragments are unique to the mutant transcripts.
#  ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)
#  

## ---- label = 'Converting gene symbols', eval = TRUE, results = 'hide'--------------------------------------------------------------------------------------------------------------
require(org.Hs.eg.db)

# Retrieve the all imported mutant transcripts.
x <- ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)

# Convert ENTREZ identifiers to gene symbols.
geneSymbols <- unique(AnnotationDbi::select(
  org.Hs.eg.db, keys = x$genomicVariants$geneID, 
  keytype = 'ENTREZID', columns = 'SYMBOL')
)

# Add the gene symbols back to the mutant transcript.
x$genomicVariants <- merge(
  x$genomicVariants, geneSymbols, 
  by.x = 'geneID', by.y = 'ENTREZID', all.x = TRUE
)

# Add the mutant transcripts (with symbols) back into the ProteoDiscography.
ProteoDiscography.hg19 <- ProteoDisco::setMutantTranscripts(
  ProteoDiscography.hg19, 
  transcripts = x$genomicVariants, 
  slotType = 'genomicVariants'
)

# We now have the SYMBOL column in our mutant transcripts DataFrame.
head(ProteoDisco::mutantTranscripts(ProteoDiscography.hg19)$genomicVariants)

## ---- label = 'Export mutant proteins to FASTA', eval = TRUE, results = 'hide'------------------------------------------------------------------------------------------------------

ProteoDisco::exportProteoDiscography(
  ProteoDiscography = ProteoDiscography.hg19, 
  outFile = 'example.fasta', 
  aggregateSamples = TRUE
)


## ---- label = 'Session Information', eval = TRUE------------------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

