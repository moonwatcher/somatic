Bootstrap the system
====================

**Build the mongodb indexes**
```
somatic rebuild
```

**Load a config file**
```
somatic populate ~/code/somatic/bootstrap/c57bl6.json

```

**Align the genes to the reference with up to 200bp flanking sequences**
```
somatic align --flanking 200
```

**Verify RSS signals around the genes**
```
somatic rss --flanking 60 --distance 9
```

**Create BWA alignment reference files**
```
somatic index
```

Load a sample
=============
```zcat A4G6U_b6_conv_fo_3_b.fastq.gz|somatic analyze --library spf_fo_b03t02```

Look at an alignment diagram
============================
```somatic sample --id M00595:21:000000000-A4G6U:1:1101:10960:6129 --format diagram```

```
R    name        S O G F  spf_fo_b03t02 | 1 | M00595:21:000000000-A4G6U:1:1101:10960:6129 | vh : J558.78.182 | dh : DSP2.x | jh : IGHJ4 | cdr3 : -20.1 2463 | productive
                          0   3   6   9   12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60  63  66  69  72  75  78  81  84  87  90  93  96  99  102 105 108 111 114 117 120 123 126 129 132 135 138 141 144 147 150 153 156 159 162 165 168 171 174 177 180 
                          AAG TTC AAG GGC AAG GCC ACA CTG ACT GTA GAC AAA TCC TCC AGC ACA GCC TAC ATG TTG CTC AGC AGC CTG ACC TCT GAG GAC TCT GCG GTC TAT TTC TGT GCA AGA GAG AGG GCC TAC TAT AGT AAC TAC GAT GCT ATG GAC TAC TGG GGT CAA GGA ACC TCA GTC ACC GTC TCC GCA G 
                          BBB BBF FFF FBB GGG GGG GGG 2FH HHH HHG HHH HCG GGG HHH HHH GHH HHH HHH HHH HFH HHH HGH HHH HHH HHG HGH HHH GHH GHH HHH GDG HHH HGG HHH HHH HGG F2F HE@ ?HE GFF ?FB 4FG BFH HHG HHH HGH HHF HHH HHH GEE GHG HHH HHH HGH HHH G00 EFE FEA EA1 1A1 F 
                          K   F   K   G   K   A   T   L   T   V   D   K   S   S   S   T   A   Y   M   L   L   S   S   L   T   S   E   D   S   A   V   Y   F   C   A   R   E   R   A   Y   Y   S   N   Y   D   A   M   D   Y   W   G   Q   G   T   S   V   T   V   S   A   
jh   IGHJ4       * -   F                                                                                                                                                                              --- T-- --- --- --- --- --- --- --- --- --- --- --- --- --- --- T-- - 
                                                                                                                                                                                                      -   Y   -   -   -   -   -   -   -   -   -   -   -   -   -   -   S   
vh   J558.78.182 * -   F  --- --- --- --- --- --- --- --T --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
dh   DSP2.x      * -   F                                                                                                                                                           -- --- --- --- --- 
v-d              * -                                                                                                                                                      --- --- - 
                                                                                                                                                                          NNN NNP P 
cdr3             * -                                                                                                                                          --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ```
```

Inspect a complete record
=========================
```somatic sample --id M00595:21:000000000-A4G6U:1:1101:10960:6129 --format json```

```json
{
    "_id": "598f619f040bf36784816e79",
    "body": {
        "FLAG": 16,
        "abundance": 1,
        "conserved cdr3 c": true,
        "conserved cdr3 w": true,
        "framed": true,
        "framed by": "82d332a8-a877-48ee-a8f6-7fbdd47f1dd0",
        "hit": [
            {
                "5 chew": {
                    "nucleotide": "AT",
                    "read frame": 2,
                    "strand": true
                },
                "AS": 92,
                "CIGAR": "129S52M5S",
                "FLAG": 16,
                "MAPQ": 60,
                "MD": "3T44T3",
                "NM": 2,
                "TLEN": 0,
                "XS": 38,
                "allele": "IGHJ4*01",
                "family": "IGHJ",
                "framed": true,
                "functionality": "F",
                "gene": "IGHJ4",
                "in frame": true,
                "picked": true,
                "primary": true,
                "query": {
                    "expected error": 0.1473087788135559,
                    "nucleotide": "TACGATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCGCAG",
                    "quality": "HHGHHHHGHHHFHHHHHHGEEGHGHHHHHHHGHHHHG00EFEFEAEA11A1F",
                    "read frame": 0,
                    "strand": false
                },
                "query end": 181,
                "query start": 129,
                "reference end": 113428565,
                "reference start": 113428513,
                "region": "jh",
                "score": 92,
                "strain": "c57bl6",
                "subject": {
                    "nucleotide": "TACTATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCAG",
                    "read frame": 0,
                    "strand": true
                },
                "subject end": 54,
                "subject id": "IGHJ4",
                "subject start": 2,
                "subject strand": false,
                "type": "gene segment",
                "uuid": "82d332a8-a877-48ee-a8f6-7fbdd47f1dd0",
                "valid": true
            },
            {
                "AS": 210,
                "CIGAR": "108M21S",
                "FLAG": 16,
                "MAPQ": 21,
                "MD": "23T84",
                "NM": 1,
                "TLEN": 0,
                "XS": 192,
                "allele": "J558.78.182*01",
                "family": "J558",
                "framed": true,
                "functionality": "F",
                "gene": "J558.78.182",
                "in frame": true,
                "picked": true,
                "primary": true,
                "query": {
                    "expected error": 0.03787568981918923,
                    "nucleotide": "AAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGTTGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                    "quality": "BBBBBFFFFFBBGGGGGGGGG2FHHHHHHGHHHHCGGGGHHHHHHGHHHHHHHHHHHHFHHHHHGHHHHHHHHHGHGHHHHGHHGHHHHHGDGHHHHGGHHHHHHHGG",
                    "read frame": 0,
                    "strand": false
                },
                "query end": 108,
                "query start": 0,
                "reference end": 115834057,
                "reference start": 115833949,
                "region": "vh",
                "score": 210,
                "strain": "c57bl6",
                "subject": {
                    "nucleotide": "AAGTTCAAGGGCAAGGCCACACTTACTGTAGACAAATCCTCCAGCACAGCCTACATGTTGCTCAGCAGCCTGACCTCTGAGGACTCTGCGGTCTATTTCTGTGCAAGA",
                    "read frame": 0,
                    "strand": true
                },
                "subject end": 305,
                "subject id": "J558.78.182",
                "subject start": 197,
                "subject strand": false,
                "type": "gene segment",
                "uuid": "11d9327c-8f82-4a78-a33c-523d02f4fb1f",
                "valid": true
            },
            {
                "3 chew": {
                    "nucleotide": "TAC",
                    "read frame": 1,
                    "strand": true
                },
                "AS": 28,
                "CIGAR": "7S14M",
                "FLAG": 16,
                "MAPQ": 17,
                "MD": "14",
                "NM": 0,
                "SA": "IGHD5-8,11,+,15S5M1S,0,0;",
                "TLEN": 0,
                "XS": 20,
                "allele": "DSP2.x*02",
                "family": "DSP2",
                "framed": true,
                "functionality": "F",
                "gene": "DSP2.x",
                "picked": true,
                "primary": true,
                "query": {
                    "expected error": 0.016409212106682676,
                    "nucleotide": "CCTACTATAGTAAC",
                    "quality": "HEGFF?FB4FGBFH",
                    "read frame": 2,
                    "strand": false
                },
                "query end": 129,
                "query start": 115,
                "reference end": 113461385,
                "reference start": 113461371,
                "region": "dh",
                "score": 28,
                "strain": "c57bl6",
                "subject": {
                    "nucleotide": "CCTACTATAGTAAC",
                    "read frame": 0,
                    "strand": true
                },
                "subject end": 14,
                "subject id": "DSP2.x",
                "subject start": 0,
                "subject strand": false,
                "type": "gene segment",
                "uuid": "8ec835aa-be07-46e7-95f5-dfbdd1d378c7",
                "valid": true
            },
            {
                "FLAG": 16,
                "n count": 5,
                "p count": 2,
                "palindrome": "NNNNNPP",
                "picked": true,
                "primary": true,
                "query": {
                    "expected error": 0.022523085031737236,
                    "nucleotide": "GAGAGGG",
                    "quality": "F2FHE@?",
                    "read frame": 0,
                    "strand": false
                },
                "query end": 115,
                "query start": 108,
                "region": "v-d",
                "subject id": "v-d",
                "uuid": "e5234732-2739-4f05-aeaf-7e1280a71b74",
                "valid": true
            },
            {
                "FLAG": 16,
                "charge": -20.099999999999998,
                "picked": true,
                "primary": true,
                "query": {
                    "expected error": 0.043196283158396444,
                    "nucleotide": "TGTGCAAGAGAGAGGGCCTACTATAGTAACTACGATGCTATGGACTACTGG",
                    "quality": "HHHHHHHGGF2FHE@?HEGFF?FB4FGBFHHHGHHHHGHHHFHHHHHHGEE",
                    "read frame": 0,
                    "strand": false
                },
                "query end": 150,
                "query start": 99,
                "region": "cdr3",
                "subject id": "cdr3",
                "uuid": "eca8151e-f0ee-4646-ba23-321d23a1e9d0",
                "valid": true,
                "weight": 2463
            }
        ],
        "id": "M00595:21:000000000-A4G6U:1:1101:10960:6129",
        "id comment": "1:N:0:ATGAG",
        "id sha1": "ca92c85e3589bf8d74d3de8fb5b01ae784cf087f",
        "in frame": true,
        "library": "spf_fo_b03t02",
        "n count": 5,
        "p count": 2,
        "palindromic": true,
        "premature": false,
        "primary": {
            "cdr3": "eca8151e-f0ee-4646-ba23-321d23a1e9d0",
            "dh": "8ec835aa-be07-46e7-95f5-dfbdd1d378c7",
            "jh": "82d332a8-a877-48ee-a8f6-7fbdd47f1dd0",
            "v-d": "e5234732-2739-4f05-aeaf-7e1280a71b74",
            "vh": "11d9327c-8f82-4a78-a33c-523d02f4fb1f"
        },
        "productive": true,
        "sequence": {
            "expected error": 0.22511439692865007,
            "nucleotide": "CTTACCTGCGGAGACGGTGACTGAGGTTCCTTGACCCCAGTAGTCCATAGCATCGTAGTTACTATAGTAGGCCCTCTCTCTTGCACAGAAATAGACCGCAGAGTCCTCAGAGGTCAGGCTGCTGAGCAACATGTAGGCTGTGCTGGAGGATTTGTCTACAGTCAGTGTGGCCTTGCCCTTGAACTT",
            "quality": "FFFFFF1A11AEAEFEFE00GHHHHGHHHHHHHGHGEEGHHHHHHFHHHGHHHHGHHHFBGF4BF?FFGEH?@EHF2FGGHHHHHHHGGHHHHGDGHHHHHGHHGHHHHGHGHHHHHHHHHGHHHHHFHHHHHHHHHHHHGHHHHHHGGGGCHHHHGHHHHHHF2GGGGGGGGGBBFFFFFBBBBB",
            "read frame": 0,
            "strand": true
        },
        "valid": true
    },
    "head": {
        "FLAG": 16,
        "abundance": 1,
        "conserved cdr3 c": true,
        "conserved cdr3 w": true,
        "framed": true,
        "framed by": "82d332a8-a877-48ee-a8f6-7fbdd47f1dd0",
        "id": "M00595:21:000000000-A4G6U:1:1101:10960:6129",
        "id sha1": "ca92c85e3589bf8d74d3de8fb5b01ae784cf087f",
        "in frame": true,
        "library": "spf_fo_b03t02",
        "n count": 5,
        "p count": 2,
        "palindromic": true,
        "premature": false,
        "productive": true,
        "valid": true
    }
}
```
