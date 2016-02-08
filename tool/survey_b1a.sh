#!/bin/zsh

SESSION_NAME="somatic_survey"
START_DIRECTORY="/Users/lg/thesis/plot"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n etc -d
new-window -t"$SESSION_NAME:1" -n one "\
"
tmux new-window -t"$SESSION_NAME:2" -n two "\
somatic -v debug survey --library 'c57bl6b02t01gfb1a'      --profile np33      --name 'gf_b1a_b02_np33';\
somatic -v debug survey --library 'c57bl6b02t01gfb1a'      --profile p33       --name 'gf_b1a_b02_p33';\
somatic -v debug survey --library 'c57bl6b03t01gfb1a'      --profile np33      --name 'gf_b1a_b03_np33';\
somatic -v debug survey --library 'c57bl6b03t01gfb1a'      --profile p33       --name 'gf_b1a_b03_p33';\
"
tmux new-window -t"$SESSION_NAME:3" -n three "\
somatic -v debug survey --library 'c57bl6b04t01gfb1a'      --profile np33      --name 'gf_b1a_b04_np33';\
somatic -v debug survey --library 'c57bl6b04t01gfb1a'      --profile p33       --name 'gf_b1a_b04_p33';\
somatic -v debug survey --library 'c57bl6b05t01gfb1a'      --profile np33      --name 'gf_b1a_b05_np33';\
somatic -v debug survey --library 'c57bl6b05t01gfb1a'      --profile p33       --name 'gf_b1a_b05_p33';\
"
tmux new-window -t"$SESSION_NAME:4" -n four "\
somatic -v debug survey --library 'c57bl6b01t01spfb1a'     --profile np33      --name 'spf_b1a_b01_np33';\
somatic -v debug survey --library 'c57bl6b01t01spfb1a'     --profile p33       --name 'spf_b1a_b01_p33';\
somatic -v debug survey --library 'c57bl6b02t01spfb1a'     --profile np33      --name 'spf_b1a_b02_np33';\
somatic -v debug survey --library 'c57bl6b02t01spfb1a'     --profile p33       --name 'spf_b1a_b02_p33';\
"
tmux new-window -t"$SESSION_NAME:5" -n five "\
somatic -v debug survey --library 'c57bl6b03t01spfb1a'     --profile np33      --name 'spf_b1a_b03_np33';\
somatic -v debug survey --library 'c57bl6b03t01spfb1a'     --profile p33       --name 'spf_b1a_b03_p33';\
somatic -v debug survey --library 'c57bl6b04t01spfb1a'     --profile np33      --name 'spf_b1a_b04_np33';\
somatic -v debug survey --library 'c57bl6b04t01spfb1a'     --profile p33       --name 'spf_b1a_b04_p33';\
"
tmux new-window -t"$SESSION_NAME:6" -n six "\
somatic -v debug survey --library 'c57bl6b05t01spfb1a'     --profile np33      --name 'spf_b1a_b05_np33';\
somatic -v debug survey --library 'c57bl6b05t01spfb1a'     --profile p33       --name 'spf_b1a_b05_p33';\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"

somatic -v debug surveyDensity \
gf_b1a_b02_np33 \
gf_b1a_b02_p33 \
gf_b1a_b03_np33 \
gf_b1a_b03_p33 \
gf_b1a_b04_np33 \
gf_b1a_b04_p33 \
gf_b1a_b05_np33 \
gf_b1a_b05_p33 \
spf_b1a_b01_np33 \
spf_b1a_b01_p33 \
spf_b1a_b02_np33 \
spf_b1a_b02_p33 \
spf_b1a_b03_np33 \
spf_b1a_b03_p33 \
spf_b1a_b04_np33 \
spf_b1a_b04_p33 \
spf_b1a_b05_np33 \
spf_b1a_b05_p33 \
> spf_gf_b1a_rep.csv && \
intensity.r spf_gf_b1a_rep.csv "SPF vs GF B1a Repetitions"

somatic -v debug surveyDensity \
gf_b1a_b02_np33 \
gf_b1a_b03_np33 \
gf_b1a_b04_np33 \
gf_b1a_b05_np33 \
spf_b1a_b01_np33 \
spf_b1a_b02_np33 \
spf_b1a_b03_np33 \
spf_b1a_b04_np33 \
spf_b1a_b05_np33 \
> spf_gf_b1a_np_rep.csv && \
intensity.r spf_gf_b1a_np_rep.csv "SPF vs GF B1a Non Productive Repetitions"

somatic -v debug surveyDensity \
spf_b1a_b01_np33 \
spf_b1a_b02_np33 \
spf_b1a_b03_np33 \
spf_b1a_b04_np33 \
spf_b1a_b05_np33 \
> gf_b1a_np_rep.csv && \
intensity.r gf_b1a_np_rep.csv "GF B1a Non Productive Repetitions"

somatic -v debug surveyDensity \
spf_b1a_b01_np33 \
spf_b1a_b02_np33 \
spf_b1a_b03_np33 \
spf_b1a_b04_np33 \
spf_b1a_b05_np33 \
> spf_b1a_np_rep.csv && \
intensity.r spf_b1a_np_rep.csv "SPF B1a Non Productive Repetitions"



somatic surveyPrevalence spf_b1a_np33 > spf_b1a_np33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_np33 > gf_b1a_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_p33 > spf_b1a_p33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_p33 > gf_b1a_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_npn133 > spf_b1a_npn133_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_npn033 > spf_b1a_npn033_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_pn133 > spf_b1a_pn133_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_pn033 > spf_b1a_pn033_vj_prevalence.csv
# somatic surveyPrevalence gf_b1a_npn133 > gf_b1a_npn133_vj_prevalence.csv
# somatic surveyPrevalence gf_b1a_npn033 > gf_b1a_npn033_vj_prevalence.csv
# somatic surveyPrevalence gf_b1a_pn133 > gf_b1a_pn133_vj_prevalence.csv
# somatic surveyPrevalence gf_b1a_pn033 > gf_b1a_pn033_vj_prevalence.csv


somatic surveyPrevalence gf_b1a_b02_np33 > gf_b1a_b02_np33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b02_p33 > gf_b1a_b02_p33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b03_np33 > gf_b1a_b03_np33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b03_p33 > gf_b1a_b03_p33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b04_np33 > gf_b1a_b04_np33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b04_p33 > gf_b1a_b04_p33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b05_np33 > gf_b1a_b05_np33_vj_prevalence.csv
somatic surveyPrevalence gf_b1a_b05_p33 > gf_b1a_b05_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b01_np33 > spf_b1a_b01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b01_p33 > spf_b1a_b01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b02_np33 > spf_b1a_b02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b02_p33 > spf_b1a_b02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b03_np33 > spf_b1a_b03_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b03_p33 > spf_b1a_b03_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b04_np33 > spf_b1a_b04_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b04_p33 > spf_b1a_b04_p33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b05_np33 > spf_b1a_b05_np33_vj_prevalence.csv
somatic surveyPrevalence spf_b1a_b05_p33 > spf_b1a_b05_p33_vj_prevalence.csv
