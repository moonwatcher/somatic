# Development of detailed high throughput analysis of V(D)J recombination

## introduction
To effectively fight off the myriad pathogens an organism encounters without exponentially increasing the size of it's genome, the mammalian adaptive immune system evolved to encode an astounding variety of antibodies from a small amount of coding sequences by rearranging the genetic makeup of each individual B cell so it can produce an antibody of different specificity. B cell antibodies are assembled from a pair of heavy chains (IgH) recombined from a V, D and J coding segments and a pair of light chains (IgL) recombined from just a V and J coding segments. The intrinsic entropy of an individual's B cell repertoire is governed by the combinatorial diversity derived from the choice of a V, D and J coding segments, nucleotide deletions and random nucleotide additions on the joining junctions and somatic hypermutations.

Recent advances in high throughput sequencing and primer design allow us to assay massive B cells populations in different species and under different conditions. However to enable rapid differential analysis our understanding of the somatic recombination process must be characterized as an equally high throughput capable software package. To that end we implemented somatic, a statistical framework for high throughput analysis of V(D)J recombinations. We used it to study the heavy chain loci of the C57BL/6 mouse strain, validate the expected behavior of the process and discover some interesting novel properties.


V, D and J coding segments are flanked by recombination signal sequences (RSS) that are recognized by the Recombination Activating Gene (RAG) proteins. An RSS is composed of a heptamer and a nonamer separated by either a 12 or a 23 base pair spacer with the heptamer side always adjacent to the coding segment. The nucleotide sequence of the heptamer and nonamer as well as the length of the spacer are essential for a successful recombination to occur and are highly conserved while the actual sequence of the spacer seems to be of lesser importance. Since somatic recombination preferentially occurs between gene segments flanked by RSSs of dissimilar length, a principal known as the 12/23 rule, the sequences recombine in the correct order.

Recombination starts with RAG nicking at the base of the heptamer between the RSS and the coding segment. Upon formation of a paired complex between the two gene segments RAG converts the nicks into double strand breaks producing two blunt ends on the side of the RSS and two hairpined ends on the side of the coding segments. The four ends are held together in a postcleavage complex by RAG. In a process though to be mediated by RAG and assisted by Artemis:DNA-PKcs complex along with Ku, DNA ligase IV and XRCC4, the two blunt signal ends are ligated directly to form a circular DNA signal joint while the hairpin ends are opened up and further processed. In the fusing process nucleotides are either added ore removed at the open ends by terminal deoxynucleotidyl transferase (TdT) in a processes that continues until there are complimentary sequences on both sides at which point and the opposite strands can pair up and and ligases can fill in the gaps.

The rearrangement process progresses step by step with with each phase being checked for success before proceeding to the next. B cell generation starts with creating pro B cell in the bone marrow in adults or in liver tissues in the fetus. During pro B cell assembly a D region will be joined by a J region on both maternal and paternal alleles. The first of the two DJ joins that forms an in frame rearrangement with a V region leads to expression of a μ heavy chain which pairs with an invariant surrogate light chain and in association with signal-transducing Igα and Igβ subunits forms the pre-B cell receptor (pre-BCR). Successful expression of a pre-BCR arrests any further attempts to assemble an IgH by down regulating RAG. Expressing a functional antibody from both maternal and paternal allele is an undesirable scenario as the secondary rearrangement could potentially produce an antibody which is self reactive but go undetected because it's signal is masked by the primary antibody.

Because of the stochastic and erroneous nature of the recombination process, only a small proportion of cells generate an in frame VDJ rearrangement, express a functional heavy chain and selected to progress to the next stage of differentiation. Expression of the pre-BCR leads to progression from pro-B cell to the pre-B cell stage where light chain gene rearrangement can take place. Pro-B cell that successfully expressed pre-BCR proliferate prior to light chain rearrangement. This insures the functional IgH rearrangement will not lost in the event of an unsuccessful light chain recombination and also further increases diversity by allowing the functional μ chain a chance to pair with multiple light chains. Cells that acquire a functional BCR, with an in frame rearranged heavy chain that successfully pairs with an in frame rearranged light chain make it out into the peripheral immune system and become mature B cells.

Because there are several members of V, D (for the heavy chain) and J segments, the combinatorial nature of V(D)J rearrangement allows generation of a vastly diverse Ig repertoire from a relatively modest amount of genetic information present in the germline. For example in ~2.7 Mb of the murine IgH locus, there are 195 VH genes, 18 DH elements and 4 JH segments (Johnston et al., 2006; Lefranc, 2003). However, there are several additional genetic mechanisms that contribute significantly to the diversity of the Ab repertoire. Among them is the imprecise processing of V, D and J segment coding ends during V(D)J recombination, and the addition of non-template encoded nucleotides (N nucleotides) at the junctions of V(D)J rearrangement. Opening of the hairpin at the end of coding joints is frequently asymmetric (Lewis, 1994), and can occur several bases away from the apex of the hairpin, resulting in overhanging nucleotides, which are then incorporated into the joint to form non-germline encoded palindromic (P) nucleotides (Lewis, 1994). In addition, during heavy chain gene rearrangement, terminal
12 deoxynucletidyl transferase (TdT) may add N nucleotides and potentially remove nucleotides from the coding sequences (Alt and Baltimore, 1982; Gilfillan et al., 1993; Komori et al., 1993; Thai et al., 2002). These processes of junctional diversity occur at the site of the Ab variable region known to make direct contact with Ag, namely the complementarity-determining region 3 (CDR3), and provide a major source of diversification for the Ig repertoire (Xu and Davis, 2000). In addition to the mechanisms that operate at the level of the Ig genes, the combination of different heavy and light chains in individual B cells further enhances the total Ab repertoire, since the variable region of both chains participate in Ag recognition. The added diversity provided by N nucleotide addition and exonuclease activity of TdT is further underscored by examination of the fetal liver B cell. In mice TdT is not expressed in fetal liver B cell progenitors (Gregoire et al., 1979; Carlsson and Holmberg, 1990) and the IgH repertoire is severely restricted since D to JH joining now occurs through the use of short sequence homologies present between the 3’ end of DH elements and 5’ end of JH elements. Since the majority of DH elements contain a stop codon in one reading frame and another reading frame results in expression of a truncated μ protein, microhomology based joining of DH and JH elements serves to insure that a single favorable DH reading frame is used. While this microhomology based joining severely restricts the IgH repertoire early in ontogeny, this process also insure a rapid production of functional V(D)J rearrangements thus establishing a first layer of immune recognition.


## Materials and Methods
## Sequence assignment
## Statistical analysis
## Shanon Entropy
## Results

## Data collection and preprocessing
To preform a high throughput analysis of a large number of samples obtained from C57BL/6 mice we started by creating an exhaustive annotated reference of the IgH loci for the strain.
The gene list provided by IMGT included annotations for functionality and reading frames for most genes but has, sadly, been very poorly annotated for strain. In an effort to assemble a coherent repertoire of the IGH genes specific to the C57BL/6 strain several other sources had to be consulted and cross referenced. There seems to be some controversy as to start location of many of the IGHV genes. Alternative sources such as MGI and ensembl seem to start many of them upstream of the start positions on IMGT. Many genes start 10 nucleotides earlier but some much more. For the purpose of this study, which focused on 146bp MiSeq read positions on the V(D)J recombination site, the actual start location of the VH genes was on little significance 
but this remains an issue to be confirmed at a later stage.

The Four immunoglobulin heavy chain joining (IGHJ) genes have been originally identified on BALB/c (Solin M.L. et al. 1992 [2]).
IGHJ4 has been identified on the C57BL/6 strain (Owens J.D. Jr et al. 1991 [3]) but the remaining 3 were not strictly annotated for C57BL/6 on any locatable accession. To established which exact variation of the remaining 3 IGHJ genes were the best match for the C57BL/6 strain the various IGHJ genes found on IMGT were used to search on accession AC073553.5, which contains the C57BL/6 IGHJ loci. IGHJ1 from accession X63164.1 from the A/J strain was an exact match to AC073553.5[133089:133142]. IGHJ2 and IGHJ3, annotated in accession X63164 for the A/J strain and accession AJ851868 for the 129S1 strain were matched perfectly to AC073553.5[133407:133455] and AC073553.5[133790:133838] respectively.

The Four immunoglobulin heavy chain diversity (IGHD) gene subgroups DFL16, DSP2, DQ52, and DST4 have been originally identified on the BALB/c mouse but previous studies have suggested significant differences between the BALB/c and C57BL/6 strains on the IGHD loci. While the BALB/c mice have at lest 13 IGHD genes, only 10 have been identified on C57BL/6, which includes one DFL16, six DSP2, one DQ52, and two DST4 genes (Jian Ye 2004 [1]. Annotated in accession AC073553.5).

The full C57BL/6 IGHV gene repertoire is not currently known. The majority of the 16 known gene families in the C57BL/6 IGHV loci has been described and annotated in accession BN000872 (Johnston C.M. et al. 2006 [4]). The IGVH genes annotated on BN000872 where than cross referenced with the gene tables on the IMGT website [5], and with the chromosome coordinates listed on MGI and ensembl. Finally the sequences from IMGT were validated against the coordinates on the BN000872 accessions and a coherent set of IGVH genes was assembled.

Several IGVH genes from alternative sources, mentioned on the IMGT gene table, mostly discovered later than the 2006 Johnston assembly, were added to that set and validated against their accessions.

Using the command line version of blat, the C57BL/6 gene set was aligned, with 200bp flanking sequences, to the GRCm38.p3 chromosome 12 reference and yielded unique and perfect alignments.

## aligning samples to the gene repertoire and picking a V/D/J match for each sample
Using the described C57BL/6 repertoire a custom database was constructed for igblast and was used to align the reads. 

## numbers
invalidate sample:
 * no igblast hits at all
 * either no valid VH leader or no valid JH leader

invalidate an igblast hit:
 * match to unknown gene (should never happen)
 * hit contains gaps
 * VH or JH hit aligns to a different strand than the framing hit
 
valid and not productive: 2511675
valid and productive: 5985117
invalid: 1228858
total: 9725650

```
 [1] Immunogenetics (2004) 56: 399–404 DOI 10.1007/s00251-004-0712-z
 [2] Solin M.L. et al., Immunogenetics, 36, 306-313 (1992). PMID: 1644448
 [3] Owens J.D. Jr et al., Mol. Cell. Biol., 11, 5660-5670 (1991). PMID: 1922069
 [4] Johnston C.M. et al., J Immunol., 176, 4221-4234 (2006). PMID: 16547259
 [5] http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=Mus_musculus&group=IGHV
```