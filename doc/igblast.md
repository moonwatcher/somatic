step 1
======

json encode the fastq read

```json
{
    "read id": "@M00595:20:000000000-A4N54:1:1101:17506:1692 1:N:0:AAGTG",
    "code": "CTTACCTGAGGAGACGGTGACTGAGGTTCCTTGACCCCAGTAGTCCATAGCATAGTAACGGATGCACAGTAATAGACCGCAGAGTCCTCAGATGTCAGGCTGCTGAGTTGCATGTAGGCTGTGCTGGAGGATTTGTCTGCAGTC",
    "quality": "BFFF?CBEFEGCBCEGEHFDGGDGHHFGHGFAGBDFFFH3GGG2FGHG5FBGFHGF5F5FHHFEF1EFFD5GF5GGGFEGE1F2FGHHGHHHH3421?BG?FGADEDGFHH4BF?CFFGEGHHFDF?EGHFGFBGHH1BBFGFH",
}
```

step 2
======

convert fastq

```
@M00595:20:000000000-A4N54:1:1101:17506:1692 1:N:0:AAGTG
CTTACCTGAGGAGACGGTGACTGAGGTTCCTTGACCCCAGTAGTCCATAGCATAGTAACGGATGCACAGTAATAGACCGCAGAGTCCTCAGATGTCAGGCTGCTGAGTTGCATGTAGGCTGTGCTGGAGGATTTGTCTGCAGTC
+
BFFF?CBEFEGCBCEGEHFDGGDGHHFGHGFAGBDFFFH3GGG2FGHG5FBGFHGF5F5FHHFEF1EFFD5GF5GGGFEGE1F2FGHHGHHHH3421?BG?FGADEDGFHH4BF?CFFGEGHHFDF?EGHFGFBGHH1BBFGFH
```

to fasta

```
>@M00595:20:000000000-A4N54:1:1101:17506:1692 1:N:0:AAGTG
AGGACTCACCTGAGGAGACTGTGAGAGTGGTGCCTTGGCCCCAGTAGTCAAAGTAGTGCTACCGTAGTAATAAGCCCCCT
CCACAGTAATAGACAGCAGGGTCCTCAGATGTCAGGCTGTTCAACACCATGTACACTGTGCTGG
```

step 3
======

execute to igblast the sequence

```
igblastn \
-germline_db_V database/mouse_gl_V \
-germline_db_J database/mouse_gl_J \
-germline_db_D database/mouse_gl_D \
-organism mouse \
-domain_system kabat \
-query query.fasta \
-auxiliary_data optional_file/mouse_gl.aux \
-show_translation \
-outfmt "7 qseqid sseqid qstart qend sstart send gapopen gaps mismatch pident bitscore evalue length sstrand"
```

return is

```
# IGBLASTN 2.2.29+
# Query: @M00595:20:000000000-A4N54:1:1101:17506:1692 1:N:0:AAGTG
# Database: database/mouse_gl_V database/mouse_gl_D database/mouse_gl_J
# Domain classification requested: kabat

# Note that your query represents the minus strand of a V gene and has been converted to the plus strand. The sequence positions refer to the converted sequence. 

# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.
J558.1.85	N/A	JH4	VH	No	Out-of-frame	No	-

# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
TGCAT	N/A	N/A	CCG	TTACT	

# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)
FR3	1	83	83	83	0	0	100
Total	N/A	N/A	83	83	0	0	100

# Hit table (the first field indicates the chain type of the hit)
# Fields: query id, subject id, q. start, q. end, s. start, s. end, gap opens, gaps, mismatches, % identity, bit score, evalue, alignment length, subject strand
# 6 hits found
V	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	J558.1.85	1	83	237	319	0	0	0	100.00	  131	7e-33	83	plus
V	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	145A	1	82	210	291	0	0	1	98.78	  126	2e-31	82	plus
V	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	H13-3	1	82	210	291	0	0	1	98.78	  126	2e-31	82	plus
J	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	JH4	87	138	2	53	0	0	0	100.00	  103	4e-27	52	plus
J	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	JH2	100	138	7	45	0	0	6	84.62	30.2	4e-05	39	plus
J	reversed|@M00595:20:000000000-A4N54:1:1101:17506:1692	JH3	103	116	13	26	0	0	1	92.86	20.3	0.043	14	plus
# BLAST processed 1 queries
```

step 4
======

convert igblast results to json and append to the node. since our match is on the minus strand
we calculate the reverse. alternatively one can switch between code and reverse and flip the 
stand attribute.

```json

{
    "read id": "@M00595:20:000000000-A4N54:1:1101:17506:1692 1:N:0:AAGTG",
    "code": "CTTACCTGAGGAGACGGTGACTGAGGTTCCTTGACCCCAGTAGTCCATAGCATAGTAACGGATGCACAGTAATAGACCGCAGAGTCCTCAGATGTCAGGCTGCTGAGTTGCATGTAGGCTGTGCTGGAGGATTTGTCTGCAGTC",
    "reverse": "GACTGCAGACAAATCCTCCAGCACAGCCTACATGCAACTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCATCCGTTACTATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCAGGTAAG",
    "quality": "BFFF?CBEFEGCBCEGEHFDGGDGHHFGHGFAGBDFFFH3GGG2FGHG5FBGFHGF5F5FHHFEF1EFFD5GF5GGGFEGE1F2FGHHGHHHH3421?BG?FGADEDGFHH4BF?CFFGEGHHFDF?EGHFGFBGHH1BBFGFH",
    "match": [
        [
            { "hit": "V,J558.1.85,1,83,237,319,0,0,0,100.00,131,7e-33,83,+,-" },
            { "hit": "V,145A,1,82,210,291,0,0,1,98.78,126,2e-31,82,+,-" },
            { "hit": "V,H13-3,1,82,210,291,0,0,1,98.78,126,2e-31,82,+,-" },
            { "hit": "J,JH4,87,138,2,53,0,0,0,100.00,103,4e-27,52,+,-" },
            { "hit": "J,JH2,100,138,7,45,0,0,6,84.62,30.2,4e-05,39,+,-" },
            { "hit": "J,JH3,103,116,13,26,0,0,1,92.86,20.3,0.043,14,+,-" }
        ]
    ]
}
```

schema is

idx |tubo                   |idx |igblast   |comment
----|-----------------------|----|----------|-------
    |query polarity         |1   |          |
    |query id               |1   |qseqid    |```(reversed|)?{igblast query id}```
0   |segment                |0   |          |
1   |subject id             |2   |sseqid    |
2   |query start            |3   |qstart    |
3   |query end              |4   |qend      |
4   |subject start          |5   |sstart    |
5   |subject end            |6   |send      |
6   |gap openings           |7   |gapopen   |
7   |gaps                   |8   |gaps      |
8   |mismatches             |9   |mismatch  |
9   |percentage identical   |10  |pident    |
10  |bit score              |11  |bitscore  |
11  |evalue                 |12  |evalue    |
12  |alignment length       |13  |length    |
13  |subject strand         |14  |sstrand   |
14  |query strand           |    |          |


regex to break an igblast hit into a dictionary

```
(?P<segment>[VDJ])\t
(?:(?P<query_polarity>reversed)\|)?(?:[^,]+)\t
(?P<subject_id>[^,]+)\t
(?P<query_start>[^,]+)\t
(?P<query_end>[^,]+)\t
(?P<subject_start>[^,]+)\t
(?P<subject_end>[^,]+)\t
(?P<gap_openings>[^,]+)\t
(?P<gaps>[^,]+)\t
(?P<mismatches>[^,]+)\t
(?P<percentage_identical>[^,]+)\t
(?P<bit_score>[^,]+)\t
(?P<evalue>[^,]+)\t
(?P<alignment_length>[^,]+)\t
(?P<subject_strand>[^,]+)
```

pattern to assemble a compressed match
```
{segment},{subject_id},{query_start},{query_end},{subject_start},{subject_end},
{gap_openings},{gaps},{mismatches},{percentage_identical},{bit_score},{evalue},
{alignment_length},{subject_strand},{query_strand}
```

a match can be expanded to the dictionary

```json
{
    "match": "V,J558.1.85,1,83,237,319,0,0,0,100.00,131,7e-33,83,+,-",
    "segment": "V",
    "subject id": "J558.8.98",
    "query start": 1,
    "query end": 83,
    "subject start": 237,
    "subject end": 319,
    "gap openings": 0,
    "gaps": 0,
    "mismatches": 0,
    "percentage identical": 100.00,
    "bit score": 131,
    "evalue": 7e-33,
    "alignment length": 83,
    "subject strand": "+",
    "query strand": "-"
}
```
