IMGT FASTA header
=================

idx |name                                               |regex
----|---------------------------------------------------|--------------------------------------------------------------------------------------
1   |IMGT accession number                              |```(P:<imgt accession>[^|]*)```
2   |gene and allele name                               |```(P:<imgt gene name>(P:<imgt subgroup>[A-Z0-9]+)-[0-9]+)\*(P:<imgt allele>[0-9]+)```
3   |species                                            |```(P:<organism name>[^|_]*)(?:_(P:<strain>[^|]+))?```
4   |functionality                                      |```(P:<functionality>F|ORF|P)```
5   |region name                                        |```(P:<region>V|D|J)-REGION```
6   |start and end positions relative to the accession  |```(P:<start>[0-9]+)\.\.(P:<start>[0-9]+)```
7   |number of nucleotides in the accession             |```(P:<length>[0-9]+)\snt```
8   |reading frame or NR for non coding                 |```(P:<reading frame>[123]|NR)```
9   |number of nucleotides added in 5'                  |```[^|]*```
10  |number of nucleotides added or removed in 3'       |```[^|]*```
11  |number of added, deleted, substituted nucleotides  |```[^|]*```
12  |number of amino acids                              |```[^|]*```
13  |number of characters in the sequence               |```[^|]*```
14  |partial                                            |```[^|]*```
15  |reverse complementary                              |```(?:rev-compl)?```

**4. functionality**

**functional** A germline entity (V-GENE, D-GENE or J-GENE) or a C-GENE is functional if the coding region has an open reading frame without stop codon, and if there is no described defect in the splicing sites, recombination signals and/or regulatory elements.

**open reading frame** A germline entity (V-GENE, D-GENE or J-GENE) or a C-GENE is qualified as ORF (Open Reading Frame) if the coding region has an open reading frame, but:
* alterations have been described in the splicing sites, recombination signals and/or regulatory elements.
* and/or changes of conserved amino acids have been suggested by the authors to lead to uncorrect folding.
* and/or the entity is an orphon.

**pseudogene** A germline entity (V-GENE, D-GENE or J-GENE) or a C-GENE is qualified as pseudogene if the coding region has stop codon(s) and/or frameshift mutation(s).
In particular, a V-GENE is considered as pseudogene if these defects occur in the L-PART1 and/or V-EXON, or if there is a mutation in the L-PART1 INIT-CODON atg.
A J-GENE is considered as pseudogene if it has been identified by the presence of a recombination signal upstream of an open reading frame, but it has no donor splicing site in 5' or the donor splice is not in the expected sf1 or if it has no conserved Phe(or Trp)-Gly-X-Gly motif. If the defects are important, pseudogenes can eventually also be qualified as vestigial (vg), for example: a germline V-GENE which cannot be assigned to a given subgroup because it is too divergent from the other pseudogenes and has too many stop codons and frameshifts.
Only vestigial pseudogenes reported in Gene tables are qualified as pseudogenes in the List of IG and TR genes and in the Potential germline repertoires.

**parentheses**, (F) and (P), when the accession number refers to a gene for which the functionality needs to be confirmed by a gDNA sequence that, for V,D,J, should be germline. For example, (F) or (P) are assigned to a C-REGION in cDNA or to a V-, D- or J-REGION in rearranged gDNA, for which the corresponding gDNA c-REGION or germline V-, D- or J-REGION has not been yet isolated.

**brackets**, [F] and [P], when the accession number refers to a gene for which the functionality needs to be confirmed by a complete sequence. For example, [F] or [P] are assigned to gDNA of V-, D- or J-REGION, not known as being germline or rearranged.


```
The FASTA header contains 15 fields separated by '|':

1. IMGT/LIGM-DB accession number(s)
2. gene and allele name
3. species
4. functionality
5. exon(s), region name(s), or extracted label(s)
6. start and end positions in the IMGT/LIGM-DB accession number(s)
7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
8. codon start, or 'NR' (not relevant) for non coding labels and out-of-frame pseudogenes
9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
12. number of amino acids (AA): this field indicates that the sequence is in amino acids
13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
14. partial (if it is)
15. reverse complementary (if it is)
```