##load##
load('somatic.gz')

This will create a *List* variable in the R session called somatic.

The List is of **S4** objects of the following deceleration:

```
setClass(
    "Survey",
    representation(
        vj="data.frame",
        vdj="data.frame",
        vj_correlation="matrix",
        vdj_correlation="list"
    )
)
```

The *vj* slot is a data frame has the following columns:

```
VJ_CDR3_charge
VJ_CDR3_length
VJ_CDR3_weight
VJ_J_5_chew
VJ_V_3_chew
VJ_V_J_N_count
VJ_V_J_P_count
VJ_V_J_length
VJ_chew
```

The *vdj* slot is a data frame has the following columns:

```
VDJ_CDR3_charge
VDJ_CDR3_length
VDJ_CDR3_weight
VDJ_D_3_chew
VDJ_D_5_chew
VDJ_D_J_N_count
VDJ_D_J_P_count
VDJ_D_J_length
VDJ_J_5_chew
VDJ_N_count
VDJ_P_count
VDJ_V_3_chew
VDJ_V_D_N_count
VDJ_V_D_P_count
VDJ_V_D_length
VDJ_chew
```

the *vj_correlation* slot is a matrix of *JH* vs *VH* correlations for the samples that had no *DH* identified.

the *vdj_correlation* slot is a List of matrices of *DH* vs *VH* correlations for each *JH*.

The matrices have the gene names as row and column labels.

