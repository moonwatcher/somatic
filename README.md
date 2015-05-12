Bootstrap the system
====================

**Build the mongodb indexes**
```
somatic rebuild
```

**Load genes for C57BL/6**
```
somatic gene-populate bootstrap/mouse_ighd.json 
somatic gene-populate bootstrap/mouse_ighj.json 
somatic gene-populate bootstrap/mouse_ighv_c57bl6_bn000872.json 
somatic gene-populate bootstrap/mouse_ighv_extra.json 
```

**Align the genes to the reference with up to 200bp flanking sequences**
```
somatic gene-align -F 200 --strain 'C57BL/6'
```

**Create FASTA files with the C57BL/6 genes**
```
somatic gene-fasta --strain 'C57BL/6' -r JH -p aligned > db/igblast/database/mouse_c57bl6_ighj
somatic gene-fasta --strain 'C57BL/6' -r VH -p aligned > db/igblast/database/mouse_c57bl6_ighv
somatic gene-fasta --strain 'C57BL/6' -r DH -p aligned > db/igblast/database/mouse_c57bl6_ighd
```

**Generate an igblast database**
```
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in db/igblast/database/mouse_c57bl6_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in db/igblast/database/mouse_c57bl6_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in db/igblast/database/mouse_c57bl6_ighd

somatic gene-igblast-aux -r VH --strain 'C57BL/6' >  db/igblast/optional_file/mouse_gl.aux
somatic gene-igblast-aux -r DH --strain 'C57BL/6' >> db/igblast/optional_file/mouse_gl.aux
somatic gene-igblast-aux -r JH --strain 'C57BL/6' >> db/igblast/optional_file/mouse_gl.aux
```

Load a sample
=============
```zcat "A4G6U l01n01 b6_spf_preb_2.fastq.gz" | somatic populate --strain C57BL/6 --library c57bl6b02t01spfpreb```

Look at an alignment diagram
============================
```somatic view -L 2 -S 128 -p productive```

```
R    name        F  c57bl6b03t03spffo | @M00595:21:000000000-A4G6U:1:1101:10492:4386 1:N:0:ATCAC | VH : J558.16.106 | DH : DSP2.2 | JH : IGHJ1 | productive | Q 36.25
                    0  2   5   8   11  14  17  20  23  26  29  32  35  38  41  44  47  50  53  56  59  62  65  68  71  74  77  80  83  86  89  92  95  98  101 104 107 110 113 116 119 122 125 128 131 134 137 140 143 
                    CC AGC ACA GCC TAC ATG GAG CTC CGC AGC CTG ACA TCT GAG GAC ACT GCA GTC TAT TAC TGT GCA AGA CGG GGT TAT GAT TAC GAC GGT TAC TGG TAC TTC GAT GTC TGG GGC ACA GGG ACC ACG GTC ACC GTC TCC TCA GGT AAG 
                    CD DFF HG< F22 F>@ 2HH G0C B>C DGF HHH HGH HGD B3F G4? F?? BCH HHH GFH HHF HHH HGE HHG GE> GGH HHH HHH HFF F>E EEA FFG HHH HHH HFG HFF HHH GGA GDH HFF GFG HHH HGH GGG HGH GGG GGG GGG GEF FCD FFF 
                       S   T   A   Y   M   E   L   R   S   L   T   S   E   D   T   A   V   Y   Y   C   A   R   R   G   Y   D   Y   D   G   Y   W   Y   F   D   V   W   G   T   G   T   T   V   T   V   S   S   G   K   
VH   J558.16.106 F  -- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
DH   DSP2.2      F                                                                                                     --- --- --- --- 
JH   IGHJ1       F                                                                                                                         --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- - 
V-D                                                                                                            --- --- 
CDR3                                                                                               --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
D-J                                                                                                                                    --- 
                                                                                                                                       PNN 
```

Inspect a complete record
=========================
```somatic info -L 1 -S 128 -p productive```
```json
{
    "_id": "55511bf4a9605a812edb5d75",
    "framed by": "IGHJ1*01-C57BL/6",
    "head": {
        "CDR3 length": 51,
        "D-J length": 3,
        "DH length": 12,
        "JH length": 52,
        "V-D length": 6,
        "VH length": 68,
        "average phred": 36.25342465753425,
        "framed": true,
        "gapped": false,
        "id": "@M00595:21:000000000-A4G6U:1:1101:10492:4386 1:N:0:ATCAC",
        "id sha1": "16b04f32de11e1ec43c5fc3d6668afd148d557d1",
        "in frame": true,
        "library": "c57bl6b03t03spffo",
        "palindromic": true,
        "premature": false,
        "productive": true,
        "region": [
            "DSP2.2*01-C57BL/6",
            "IGHJ1*01-C57BL/6",
            "D-J",
            "J558.16.106*01-C57BL/6",
            "CDR3",
            "V-D"
        ],
        "valid": true
    },
    "hit": [
        {
            "alignment length": 68,
            "allele": "J558.16.106*01",
            "average phread": 34.64705882352941,
            "bit score": 107.0,
            "evalue": 3e-26,
            "family": "J558",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "J558.16.106",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": true,
            "query": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "quality": "CDDFFHG<F22F>@2HHG0CB>CDGFHHHHGHHGDB3FG4?F??BCHHHHGFHHHFHHHHGEHHGGE>",
                "read frame": 2,
                "strand": false
            },
            "query end": 68,
            "query start": 0,
            "region": "VH",
            "score": 107.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": true
            },
            "subject end": 294,
            "subject id": "J558.16.106*01-C57BL/6",
            "subject start": 226,
            "subject strand": true,
            "uuid": "010900c4-d97a-4c73-96f8-cbd9ec702749",
            "valid": true
        },
        {
            "alignment length": 68,
            "allele": "J558.26.116*01",
            "bit score": 104.0,
            "evalue": 2e-25,
            "family": "J558",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "J558.26.116",
            "identical": 0.9853000000000001,
            "in frame": true,
            "mismatch": 1,
            "picked": false,
            "query": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "quality": "CDDFFHG<F22F>@2HHG0CB>CDGFHHHHGHHGDB3FG4?F??BCHHHHGFHHHFHHHHGEHHGGE>",
                "read frame": 2,
                "strand": false
            },
            "query end": 68,
            "query start": 0,
            "region": "VH",
            "score": 104.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": true
            },
            "subject end": 294,
            "subject id": "J558.26.116*01-C57BL/6",
            "subject start": 226,
            "subject strand": true,
            "uuid": "34f6deb0-a6bb-4693-985a-8185bcd7bdf9",
            "valid": true
        },
        {
            "alignment length": 68,
            "allele": "J558.34.124*01",
            "bit score": 104.0,
            "evalue": 2e-25,
            "family": "J558",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "J558.34.124",
            "identical": 0.9853000000000001,
            "in frame": true,
            "mismatch": 1,
            "picked": false,
            "query": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "quality": "CDDFFHG<F22F>@2HHG0CB>CDGFHHHHGHHGDB3FG4?F??BCHHHHGFHHHFHHHHGEHHGGE>",
                "read frame": 2,
                "strand": false
            },
            "query end": 68,
            "query start": 0,
            "region": "VH",
            "score": 104.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": true
            },
            "subject end": 294,
            "subject id": "J558.34.124*01-C57BL/6",
            "subject start": 226,
            "subject strand": true,
            "uuid": "221667e0-93bc-4ca4-af52-8f63510fff0d",
            "valid": true
        },
        {
            "5 chew": {
                "nucleotide": "TCTAC",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 12,
            "allele": "DSP2.2*01",
            "average phread": 36.333333333333336,
            "bit score": 24.4,
            "evalue": 0.001,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.2",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 0,
            "picked": true,
            "query": {
                "nucleotide": "TATGATTACGAC",
                "quality": "HHHHFFF>EEEA",
                "read frame": 0,
                "strand": false
            },
            "query end": 86,
            "query start": 74,
            "region": "DH",
            "score": 24.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TATGATTACGAC",
                "read frame": 1,
                "strand": true
            },
            "subject end": 17,
            "subject id": "DSP2.2*01-C57BL/6",
            "subject start": 5,
            "subject strand": true,
            "uuid": "5a5e18ec-047b-4f58-ba86-656d87d5cbc1",
            "valid": true
        },
        {
            "3 chew": {
                "nucleotide": "TACTAC",
                "read frame": 1,
                "strand": true
            },
            "5 chew": {
                "nucleotide": "TCTATGAT",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 3,
            "allele": "DSP2.9*01",
            "bit score": 14.4,
            "evalue": 1.3,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.9",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 4,
            "picked": false,
            "query": {
                "nucleotide": "GGT",
                "quality": "FFG",
                "read frame": 0,
                "strand": false
            },
            "query end": 89,
            "query start": 86,
            "region": "DH",
            "score": 12.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "GGT",
                "read frame": 1,
                "strand": true
            },
            "subject end": 11,
            "subject id": "DSP2.9*01-C57BL/6",
            "subject start": 8,
            "subject strand": true,
            "uuid": "b2f22881-0fa5-4927-bada-b30d3c04f676",
            "valid": true
        },
        {
            "5 chew": {
                "nucleotide": "TCTACTATGG",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 7,
            "allele": "DSP2.3*01",
            "bit score": 14.4,
            "evalue": 1.3,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.3",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 0,
            "picked": false,
            "query": {
                "nucleotide": "TTACGAC",
                "quality": "FF>EEEA",
                "read frame": 1,
                "strand": false
            },
            "query end": 86,
            "query start": 79,
            "region": "DH",
            "score": 14.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TTACGAC",
                "read frame": 2,
                "strand": true
            },
            "subject end": 17,
            "subject id": "DSP2.3*01-C57BL/6",
            "subject start": 10,
            "subject strand": true,
            "uuid": "d0147038-2bf9-461f-9b68-80780c9c9053",
            "valid": true
        },
        {
            "3 chew": {
                "nucleotide": "GGTTACTAC",
                "read frame": 1,
                "strand": true
            },
            "5 chew": {
                "nucleotide": "TC",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 6,
            "allele": "DSP2.9*01",
            "bit score": 12.4,
            "evalue": 5.2,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.9",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 0,
            "picked": false,
            "query": {
                "nucleotide": "TATGAT",
                "quality": "HHHHFF",
                "read frame": 0,
                "strand": false
            },
            "query end": 80,
            "query start": 74,
            "region": "DH",
            "score": 12.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TATGAT",
                "read frame": 1,
                "strand": true
            },
            "subject end": 8,
            "subject id": "DSP2.9*01-C57BL/6",
            "subject start": 2,
            "subject strand": true,
            "uuid": "dbe155c4-1c65-4f3b-ad3b-5439840840fb",
            "valid": true
        },
        {
            "3 chew": {
                "nucleotide": "TACGAC",
                "read frame": 1,
                "strand": true
            },
            "5 chew": {
                "nucleotide": "TCTACTAT",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 3,
            "allele": "DSP2.3*01",
            "bit score": 12.4,
            "evalue": 5.2,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.3",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 3,
            "picked": false,
            "query": {
                "nucleotide": "GGT",
                "quality": "FFG",
                "read frame": 0,
                "strand": false
            },
            "query end": 89,
            "query start": 86,
            "region": "DH",
            "score": 10.9,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "GGT",
                "read frame": 1,
                "strand": true
            },
            "subject end": 11,
            "subject id": "DSP2.3*01-C57BL/6",
            "subject start": 8,
            "subject strand": true,
            "uuid": "f097bf08-a2c4-4e77-8ade-6fb3d40f78d3",
            "valid": true
        },
        {
            "5 chew": {
                "nucleotide": "C",
                "read frame": 1,
                "strand": true
            },
            "alignment length": 52,
            "allele": "IGHJ1*01",
            "average phread": 38.01923076923077,
            "bit score": 103.0,
            "evalue": 2e-27,
            "family": "IGHJ",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "IGHJ1",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": true,
            "query": {
                "nucleotide": "TACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG",
                "quality": "HHHHHHHFGHFFHHHGGAGDHHFFGFGHHHHGHGGGHGHGGGGGGGGGGEFF",
                "read frame": 0,
                "strand": false
            },
            "query end": 141,
            "query start": 89,
            "region": "JH",
            "score": 103.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG",
                "read frame": 0,
                "strand": true
            },
            "subject end": 53,
            "subject id": "IGHJ1*01-C57BL/6",
            "subject start": 1,
            "subject strand": true,
            "uuid": "d1fe985e-b81a-4d58-a60d-daa858a44969",
            "valid": true
        },
        {
            "5 chew": {
                "nucleotide": "CCTGGTTTGCTTA",
                "read frame": 2,
                "strand": true
            },
            "alignment length": 7,
            "allele": "IGHJ3*01",
            "bit score": 14.4,
            "evalue": 1.1,
            "family": "IGHJ",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "IGHJ3",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": false,
            "query": {
                "nucleotide": "CTGGGGC",
                "quality": "AGDHHFF",
                "read frame": 1,
                "strand": false
            },
            "query end": 113,
            "query start": 106,
            "region": "JH",
            "score": 14.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "CTGGGGC",
                "read frame": 1,
                "strand": true
            },
            "subject end": 20,
            "subject id": "IGHJ3*01-C57BL/6",
            "subject start": 13,
            "subject strand": true,
            "uuid": "83c1b70a-177c-4c64-a2ef-5de07b61c4e5",
            "valid": true
        },
        {
            "5 chew": {
                "nucleotide": "CCTGGTTTGC",
                "read frame": 2,
                "strand": true
            },
            "alignment length": 7,
            "allele": "IGHJ3*01",
            "bit score": 14.4,
            "evalue": 1.1,
            "family": "IGHJ",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "IGHJ3",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": false,
            "query": {
                "nucleotide": "TTACTGG",
                "quality": "GHHHHHH",
                "read frame": 1,
                "strand": false
            },
            "query end": 95,
            "query start": 88,
            "region": "JH",
            "score": 14.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TTACTGG",
                "read frame": 1,
                "strand": true
            },
            "subject end": 17,
            "subject id": "IGHJ3*01-C57BL/6",
            "subject start": 10,
            "subject strand": true,
            "uuid": "dc9be552-684d-452b-a703-be59edb48d97",
            "valid": true
        }
    ],
    "matched": true,
    "region": {
        "CDR3": {
            "average phread": 37.450980392156865,
            "picked": true,
            "query": {
                "nucleotide": "TGTGCAAGACGGGGTTATGATTACGACGGTTACTGGTACTTCGATGTCTGG",
                "quality": "HGEHHGGE>GGHHHHHHHHFFF>EEEAFFGHHHHHHHFGHFFHHHGGAGDH",
                "read frame": 0,
                "strand": false
            },
            "query end": 110,
            "query start": 59,
            "region": "CDR3",
            "subject id": "CDR3",
            "subject strand": true
        },
        "D-J": {
            "average phread": 37.333333333333336,
            "palindrome": "PNN",
            "palindrome ratio": 2.0,
            "picked": true,
            "query": {
                "nucleotide": "GGT",
                "quality": "FFG",
                "read frame": 0,
                "strand": false
            },
            "query end": 89,
            "query start": 86,
            "region": "D-J",
            "subject id": "D-J",
            "subject strand": true
        },
        "DH": {
            "5 chew": {
                "nucleotide": "TCTAC",
                "read frame": 0,
                "strand": true
            },
            "alignment length": 12,
            "allele": "DSP2.2*01",
            "average phread": 36.333333333333336,
            "bit score": 24.4,
            "evalue": 0.001,
            "family": "DSP2",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "DSP2.2",
            "identical": 1.0,
            "in frame": false,
            "mismatch": 0,
            "overlap": 0,
            "picked": true,
            "query": {
                "nucleotide": "TATGATTACGAC",
                "quality": "HHHHFFF>EEEA",
                "read frame": 0,
                "strand": false
            },
            "query end": 86,
            "query start": 74,
            "region": "DH",
            "score": 24.4,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TATGATTACGAC",
                "read frame": 1,
                "strand": true
            },
            "subject end": 17,
            "subject id": "DSP2.2*01-C57BL/6",
            "subject start": 5,
            "subject strand": true,
            "uuid": "5a5e18ec-047b-4f58-ba86-656d87d5cbc1",
            "valid": true
        },
        "JH": {
            "5 chew": {
                "nucleotide": "C",
                "read frame": 1,
                "strand": true
            },
            "alignment length": 52,
            "allele": "IGHJ1*01",
            "average phread": 38.01923076923077,
            "bit score": 103.0,
            "evalue": 2e-27,
            "family": "IGHJ",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "IGHJ1",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": true,
            "query": {
                "nucleotide": "TACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG",
                "quality": "HHHHHHHFGHFFHHHGGAGDHHFFGFGHHHHGHGGGHGHGGGGGGGGGGEFF",
                "read frame": 0,
                "strand": false
            },
            "query end": 141,
            "query start": 89,
            "region": "JH",
            "score": 103.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "TACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG",
                "read frame": 0,
                "strand": true
            },
            "subject end": 53,
            "subject id": "IGHJ1*01-C57BL/6",
            "subject start": 1,
            "subject strand": true,
            "uuid": "d1fe985e-b81a-4d58-a60d-daa858a44969",
            "valid": true
        },
        "V-D": {
            "average phread": 38.666666666666664,
            "palindrome": "NNNNNN",
            "palindrome ratio": 1.0,
            "picked": true,
            "query": {
                "nucleotide": "CGGGGT",
                "quality": "GGHHHH",
                "read frame": 0,
                "strand": false
            },
            "query end": 74,
            "query start": 68,
            "region": "V-D",
            "subject id": "V-D",
            "subject strand": true
        },
        "VH": {
            "alignment length": 68,
            "allele": "J558.16.106*01",
            "average phread": 34.64705882352941,
            "bit score": 107.0,
            "evalue": 3e-26,
            "family": "J558",
            "framed": true,
            "functionality": "F",
            "gap openings": 0,
            "gapped": false,
            "gaps": 0,
            "gene": "J558.16.106",
            "identical": 1.0,
            "in frame": true,
            "mismatch": 0,
            "picked": true,
            "query": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "quality": "CDDFFHG<F22F>@2HHG0CB>CDGFHHHHGHHGDB3FG4?F??BCHHHHGFHHHFHHHHGEHHGGE>",
                "read frame": 2,
                "strand": false
            },
            "query end": 68,
            "query start": 0,
            "region": "VH",
            "score": 107.0,
            "strain": "C57BL/6",
            "subject": {
                "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGA",
                "read frame": 2,
                "strand": true
            },
            "subject end": 294,
            "subject id": "J558.16.106*01-C57BL/6",
            "subject start": 226,
            "subject strand": true,
            "uuid": "010900c4-d97a-4c73-96f8-cbd9ec702749",
            "valid": true
        }
    },
    "sequence": {
        "nucleotide": "CCAGCACAGCCTACATGGAGCTCCGCAGCCTGACATCTGAGGACACTGCAGTCTATTACTGTGCAAGACGGGGTTATGATTACGACGGTTACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAGGTAAG",
        "quality": "CDDFFHG<F22F>@2HHG0CB>CDGFHHHHGHHGDB3FG4?F??BCHHHHGFHHHFHHHHGEHHGGE>GGHHHHHHHHFFF>EEEAFFGHHHHHHHFGHFFHHHGGAGDHHFFGFGHHHHGHGGGHGHGGGGGGGGGGEFFCDFFF",
        "read frame": 2,
        "strand": false
    },
    "strand": true
}
```
