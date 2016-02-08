#!/bin/bash

somatic -v debug surveyDensity spf_np33            gf_np33         spf_p33         gf_p33          > spf_gf.csv            && intensity.r spf_gf.csv            "SPF vs GF"
somatic -v debug surveyDensity spf_b1a_np33        gf_b1a_np33     spf_b1a_p33     gf_b1a_p33      > spf_gf_b1a.csv        && intensity.r spf_gf_b1a.csv        "B1a SPF vs GF"
somatic -v debug surveyDensity spf_fo_np33         gf_fo_np33      spf_fo_p33      gf_fo_p33       > spf_gf_fo.csv         && intensity.r spf_gf_fo.csv         "Follicular spleen SPF vs GF"
somatic -v debug surveyDensity spf_mz_np33         gf_mz_np33      spf_mz_p33      gf_mz_p33       > spf_gf_mz.csv         && intensity.r spf_gf_mz.csv         "MZ SPF vs GF"

somatic -v debug surveyDensity spf_b1a_np33        spf_iftl_np33   spf_b1a_p33     spf_iftl_p33    > spf_b1a_iftl.csv      && intensity.r spf_b1a_iftl.csv      "SPF B1a vs Immature Fetal Liver"
somatic -v debug surveyDensity spf_b1a_np33        spf_fo_np33     spf_b1a_p33     spf_fo_p33      > spf_b1a_fo.csv        && intensity.r spf_b1a_fo.csv        "SPF B1a vs Follicular spleen"
somatic -v debug surveyDensity spf_b1a_npn133      spf_b1a_npn033  spf_b1a_pn133   spf_b1a_pn033   > spf_b1a_n0_n1.csv     && intensity.r spf_b1a_n0_n1.csv     "SPF B1a N0 vs N1"
somatic -v debug surveyDensity spf_b1a_npn033      spf_iftl_np33   spf_b1a_pn033   spf_iftl_p33    > spf_b1a_n0_iftl.csv   && intensity.r spf_b1a_n0_iftl.csv   "SPF B1a N0 vs Immature Fetal Liver"
somatic -v debug surveyDensity spf_b1a_npn133      spf_iftl_np33   spf_b1a_pn133   spf_iftl_p33    > spf_b1a_n1_iftl.csv   && intensity.r spf_b1a_n1_iftl.csv   "SPF B1a N1 vs Immature Fetal Liver"

somatic -v debug surveyDensity spf_fo_np33         spf_b1a_npn033  spf_fo_p33      spf_b1a_pn033   > spf_fo_b1a_n0.csv     && intensity.r spf_fo_b1a_n0.csv     "SPF B1a N0 vs Follicular spleen"
somatic -v debug surveyDensity spf_fo_np33         spf_iftl_np33   spf_fo_p33      spf_iftl_p33    > spf_fo_iftl.csv       && intensity.r spf_fo_iftl.csv       "SPF Follicular spleen vs Immature Fetal Liver"

somatic -v debug surveyDensity spf_preb_np33       spf_b1a_np33    spf_preb_p33    spf_b1a_p33     > spf_preb_b1a.csv      && intensity.r spf_preb_b1a.csv      "SPF PreB vs B1a"
somatic -v debug surveyDensity spf_preb_np33       spf_fo_np33     spf_preb_p33    spf_fo_p33      > spf_preb_fo.csv       && intensity.r spf_preb_fo.csv       "SPF PreB vs Follicular spleen"
somatic -v debug surveyDensity spf_preb_np33       spf_mz_np33     spf_preb_p33    spf_mz_p33      > spf_preb_mz.csv       && intensity.r spf_preb_mz.csv       "SPF PreB vs MZ"

somatic -v debug surveyDensity spf_prebabm_np33    spf_prebabm_p33                                 > spf_prebabm.csv       && intensity.r spf_prebabm.csv       "SPF PreB Adult Bone Marrow"
somatic -v debug surveyDensity spf_prebftl_np33    spf_prebftl_p33                                 > spf_prebftl.csv       && intensity.r spf_prebftl.csv       "SPF PreB Fetal Liver"
somatic -v debug surveyDensity spf_iabm_np33       spf_iabm_p33                                    > spf_iabm.csv          && intensity.r spf_iabm.csv          "SPF Immature Adult Bone Marrow"
somatic -v debug surveyDensity spf_iaspl_np33      spf_iaspl_p33                                   > spf_iaspl.csv         && intensity.r spf_iaspl.csv         "SPF Immature Adult Spleen"
somatic -v debug surveyDensity spf_iftl_np33       spf_iftl_p33                                    > spf_iftl.csv          && intensity.r spf_iftl.csv          "SPF Immature Fetal Liver"


somatic -v debug surveyDensity spf_preb_np33       spf_b1a_npn033  spf_preb_p33    spf_b1a_pn033> spf_preb_b1a_n0.csv     && intensity.r spf_preb_b1a_n0.csv     "SPF PreB vs B1a N0"
somatic -v debug surveyDensity spf_preb_np33       spf_b1a_npn133  spf_preb_p33    spf_b1a_pn133> spf_preb_b1a_n1.csv     && intensity.r spf_preb_b1a_n1.csv     "SPF PreB vs B1a N1"

somatic -v debug surveyDensity spf_preb_np33       gf_preb_np33    spf_preb_p33    gf_preb_p33> spf_gf_preb.csv           && intensity.r spf_gf_preb.csv     "PreB SPF vs GF"




