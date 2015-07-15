# Development of detailed high throughput analysis of V(D)J recombination


## introduction
## Materials and Methods
## Sequence assignment
## Statistical analysis
## Shanon Entropy
## Results


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

## assembling a heavy chain gene repertoire for the C57BL/6 mouse strain

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


```
 [1] Immunogenetics (2004) 56: 399–404 DOI 10.1007/s00251-004-0712-z
 [2] Solin M.L. et al., Immunogenetics, 36, 306-313 (1992). PMID: 1644448
 [3] Owens J.D. Jr et al., Mol. Cell. Biol., 11, 5660-5670 (1991). PMID: 1922069
 [4] Johnston C.M. et al., J Immunol., 176, 4221-4234 (2006). PMID: 16547259
 [5] http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=Mus_musculus&group=IGHV
```