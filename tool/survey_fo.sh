#!/bin/zsh

SESSION_NAME="somatic_survey"
START_DIRECTORY="/Users/lg/thesis/plot"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n etc -d
new-window -t"$SESSION_NAME:1" -n one "\
"
tmux new-window -t"$SESSION_NAME:2" -n two "\
somatic -v debug survey --library 'c57bl6b02t02spffo'     --profile p33       --name 'spf_fo_b02_t02_p33';\
somatic -v debug survey --library 'c57bl6b03t01spffo'     --profile np33      --name 'spf_fo_b03_t01_np33';\
somatic -v debug survey --library 'c57bl6b03t01spffo'     --profile p33       --name 'spf_fo_b03_t01_p33';\
somatic -v debug survey --library 'c57bl6b03t02spffo'     --profile np33      --name 'spf_fo_b03_t02_np33';\
somatic -v debug survey --library 'c57bl6b03t02spffo'     --profile p33       --name 'spf_fo_b03_t02_p33';\
somatic -v debug survey --library 'c57bl6b03t03spffo'     --profile np33      --name 'spf_fo_b03_t03_np33';\
somatic -v debug survey --library 'c57bl6b03t03spffo'     --profile p33       --name 'spf_fo_b03_t03_p33';\
somatic -v debug survey --library 'c57bl6b04t01spffo'     --profile np33      --name 'spf_fo_b04_t01_np33';\
somatic -v debug survey --library 'c57bl6b04t01spffo'     --profile p33       --name 'spf_fo_b04_t01_p33';\
"
tmux new-window -t"$SESSION_NAME:3" -n three "\
somatic -v debug survey --library 'c57bl6b04t02spffo'     --profile np33      --name 'spf_fo_b04_t02_np33';\
somatic -v debug survey --library 'c57bl6b04t02spffo'     --profile p33       --name 'spf_fo_b04_t02_p33';\
somatic -v debug survey --library 'c57bl6b05t01spffo'     --profile np33      --name 'spf_fo_b05_t01_np33';\
somatic -v debug survey --library 'c57bl6b05t01spffo'     --profile p33       --name 'spf_fo_b05_t01_p33';\
somatic -v debug survey --library 'c57bl6b05t02spffo'     --profile np33      --name 'spf_fo_b05_t02_np33';\
somatic -v debug survey --library 'c57bl6b05t02spffo'     --profile p33       --name 'spf_fo_b05_t02_p33';\
somatic -v debug survey --library 'c57bl6b05t03spffo'     --profile np33      --name 'spf_fo_b05_t03_np33';\
somatic -v debug survey --library 'c57bl6b05t03spffo'     --profile p33       --name 'spf_fo_b05_t03_p33';\
somatic -v debug survey --library 'c57bl6b05t04spffo'     --profile np33      --name 'spf_fo_b05_t04_np33';\
"
tmux new-window -t"$SESSION_NAME:4" -n four "\
somatic -v debug survey --library 'c57bl6b05t04spffo'     --profile p33       --name 'spf_fo_b05_t04_p33';\
somatic -v debug survey --library 'c57bl6b02t01gffo'     --profile np33      --name 'gf_fo_b02_t01_np33';\
somatic -v debug survey --library 'c57bl6b02t01gffo'     --profile p33       --name 'gf_fo_b02_t01_p33';\
somatic -v debug survey --library 'c57bl6b02t02gffo'     --profile np33      --name 'gf_fo_b02_t02_np33';\
somatic -v debug survey --library 'c57bl6b02t02gffo'     --profile p33       --name 'gf_fo_b02_t02_p33';\
somatic -v debug survey --library 'c57bl6b02t03gffo'     --profile np33      --name 'gf_fo_b02_t03_np33';\
somatic -v debug survey --library 'c57bl6b02t03gffo'     --profile p33       --name 'gf_fo_b02_t03_p33';\
somatic -v debug survey --library 'c57bl6b03t01gffo'     --profile np33      --name 'gf_fo_b03_t01_np33';\
somatic -v debug survey --library 'c57bl6b03t01gffo'     --profile p33       --name 'gf_fo_b03_t01_p33';\
"
tmux new-window -t"$SESSION_NAME:5" -n five "\
somatic -v debug survey --library 'c57bl6b03t02gffo'     --profile np33      --name 'gf_fo_b03_t02_np33';\
somatic -v debug survey --library 'c57bl6b03t02gffo'     --profile p33       --name 'gf_fo_b03_t02_p33';\
somatic -v debug survey --library 'c57bl6b03t03gffo'     --profile np33      --name 'gf_fo_b03_t03_np33';\
somatic -v debug survey --library 'c57bl6b03t03gffo'     --profile p33       --name 'gf_fo_b03_t03_p33';\
somatic -v debug survey --library 'c57bl6b03t04gffo'     --profile np33      --name 'gf_fo_b03_t04_np33';\
somatic -v debug survey --library 'c57bl6b03t04gffo'     --profile p33       --name 'gf_fo_b03_t04_p33';\
somatic -v debug survey --library 'c57bl6b04t01gffo'     --profile np33      --name 'gf_fo_b04_t01_np33';\
somatic -v debug survey --library 'c57bl6b04t01gffo'     --profile p33       --name 'gf_fo_b04_t01_p33';\
somatic -v debug survey --library 'c57bl6b04t02gffo'     --profile np33      --name 'gf_fo_b04_t02_np33';\
"
tmux new-window -t"$SESSION_NAME:6" -n six "\
somatic -v debug survey --library 'c57bl6b04t02gffo'     --profile p33       --name 'gf_fo_b04_t02_p33';\
somatic -v debug survey --library 'c57bl6b04t03gffo'     --profile np33      --name 'gf_fo_b04_t03_np33';\
somatic -v debug survey --library 'c57bl6b04t03gffo'     --profile p33       --name 'gf_fo_b04_t03_p33';\
somatic -v debug survey --library 'c57bl6b04t04gffo'     --profile np33      --name 'gf_fo_b04_t04_np33';\
somatic -v debug survey --library 'c57bl6b04t04gffo'     --profile p33       --name 'gf_fo_b04_t04_p33';\
somatic -v debug survey --library 'c57bl6b05t01gffo'     --profile np33      --name 'gf_fo_b05_t01_np33';\
somatic -v debug survey --library 'c57bl6b05t01gffo'     --profile p33       --name 'gf_fo_b05_t01_p33';\
somatic -v debug survey --library 'c57bl6b05t02gffo'     --profile np33      --name 'gf_fo_b05_t02_np33';\
somatic -v debug survey --library 'c57bl6b05t02gffo'     --profile p33       --name 'gf_fo_b05_t02_p33';\
somatic -v debug survey --library 'c57bl6b05t03gffo'     --profile np33      --name 'gf_fo_b05_t03_np33';\
somatic -v debug survey --library 'c57bl6b05t03gffo'     --profile p33       --name 'gf_fo_b05_t03_p33';\
"
tmux new-window -t"$SESSION_NAME:7" -n seven "\
somatic -v debug survey --library 'c57bl6b01t01spffo'     --profile np33      --name 'spf_fo_b01_t01_np33';\
somatic -v debug survey --library 'c57bl6b01t01spffo'     --profile p33       --name 'spf_fo_b01_t01_p33';\
somatic -v debug survey --library 'c57bl6b01t02spffo'     --profile np33      --name 'spf_fo_b01_t02_np33';\
somatic -v debug survey --library 'c57bl6b01t02spffo'     --profile p33       --name 'spf_fo_b01_t02_p33';\
somatic -v debug survey --library 'c57bl6b01t03spffo'     --profile np33      --name 'spf_fo_b01_t03_np33';\
somatic -v debug survey --library 'c57bl6b01t03spffo'     --profile p33       --name 'spf_fo_b01_t03_p33';\
somatic -v debug survey --library 'c57bl6b02t01spffo'     --profile np33      --name 'spf_fo_b02_t01_np33';\
somatic -v debug survey --library 'c57bl6b02t01spffo'     --profile p33       --name 'spf_fo_b02_t01_p33';\
somatic -v debug survey --library 'c57bl6b02t02spffo'     --profile np33      --name 'spf_fo_b02_t02_np33';\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"

somatic -v debug surveyDensity \
spf_fo_b01_t01_np33 \
spf_fo_b01_t01_p33 \
spf_fo_b01_t02_np33 \
spf_fo_b01_t02_p33 \
spf_fo_b01_t03_np33 \
spf_fo_b01_t03_p33 \
spf_fo_b02_t01_np33 \
spf_fo_b02_t01_p33 \
spf_fo_b02_t02_np33 \
spf_fo_b02_t02_p33 \
spf_fo_b03_t01_np33 \
spf_fo_b03_t01_p33 \
spf_fo_b03_t02_np33 \
spf_fo_b03_t02_p33 \
spf_fo_b03_t03_np33 \
spf_fo_b03_t03_p33 \
spf_fo_b04_t01_np33 \
spf_fo_b04_t01_p33 \
spf_fo_b04_t02_np33 \
spf_fo_b04_t02_p33 \
spf_fo_b05_t01_np33 \
spf_fo_b05_t01_p33 \
spf_fo_b05_t02_np33 \
spf_fo_b05_t02_p33 \
spf_fo_b05_t03_np33 \
spf_fo_b05_t03_p33 \
spf_fo_b05_t04_np33 \
spf_fo_b05_t04_p33 \
gf_fo_b02_t01_np33 \
gf_fo_b02_t01_p33 \
gf_fo_b02_t02_np33 \
gf_fo_b02_t02_p33 \
gf_fo_b02_t03_np33 \
gf_fo_b02_t03_p33 \
gf_fo_b03_t01_np33 \
gf_fo_b03_t01_p33 \
gf_fo_b03_t02_np33 \
gf_fo_b03_t02_p33 \
gf_fo_b03_t03_np33 \
gf_fo_b03_t03_p33 \
gf_fo_b03_t04_np33 \
gf_fo_b03_t04_p33 \
gf_fo_b04_t01_np33 \
gf_fo_b04_t01_p33 \
gf_fo_b04_t02_np33 \
gf_fo_b04_t02_p33 \
gf_fo_b04_t03_np33 \
gf_fo_b04_t03_p33 \
gf_fo_b04_t04_np33 \
gf_fo_b04_t04_p33 \
gf_fo_b05_t01_np33 \
gf_fo_b05_t01_p33 \
gf_fo_b05_t02_np33 \
gf_fo_b05_t02_p33 \
gf_fo_b05_t03_np33 \
gf_fo_b05_t03_p33 \
> spf_gf_fo_rep.csv && \
intensity.r spf_gf_fo_rep.csv "SPF vs GF FO Repetitions"

somatic -v debug surveyDensity \
spf_fo_b01_t01_np33 \
spf_fo_b01_t02_np33 \
spf_fo_b01_t03_np33 \
spf_fo_b02_t01_np33 \
spf_fo_b02_t02_np33 \
spf_fo_b03_t01_np33 \
spf_fo_b03_t02_np33 \
spf_fo_b03_t03_np33 \
spf_fo_b04_t01_np33 \
spf_fo_b04_t02_np33 \
spf_fo_b05_t01_np33 \
spf_fo_b05_t02_np33 \
spf_fo_b05_t03_np33 \
spf_fo_b05_t04_np33 \
gf_fo_b02_t01_np33 \
gf_fo_b02_t02_np33 \
gf_fo_b02_t03_np33 \
gf_fo_b03_t01_np33 \
gf_fo_b03_t02_np33 \
gf_fo_b03_t03_np33 \
gf_fo_b03_t04_np33 \
gf_fo_b04_t01_np33 \
gf_fo_b04_t02_np33 \
gf_fo_b04_t03_np33 \
gf_fo_b04_t04_np33 \
gf_fo_b05_t01_np33 \
gf_fo_b05_t02_np33 \
gf_fo_b05_t03_np33 \
> spf_gf_fo_np_rep.csv && \
intensity.r spf_gf_fo_np_rep.csv "SPF vs GF FO Non Productive Repetitions"

somatic -v debug surveyDensity \
gf_fo_b02_t01_np33 \
gf_fo_b02_t02_np33 \
gf_fo_b02_t03_np33 \
gf_fo_b03_t01_np33 \
gf_fo_b03_t02_np33 \
gf_fo_b03_t03_np33 \
gf_fo_b03_t04_np33 \
gf_fo_b04_t01_np33 \
gf_fo_b04_t02_np33 \
gf_fo_b04_t03_np33 \
gf_fo_b04_t04_np33 \
gf_fo_b05_t01_np33 \
gf_fo_b05_t02_np33 \
gf_fo_b05_t03_np33 \
> gf_fo_np_rep.csv && \
intensity.r gf_fo_np_rep.csv "SPF vs GF FO Non Productive Repetitions"

somatic -v debug surveyDensity \
spf_fo_b01_t01_np33 \
spf_fo_b01_t02_np33 \
spf_fo_b01_t03_np33 \
spf_fo_b02_t01_np33 \
spf_fo_b02_t02_np33 \
spf_fo_b03_t01_np33 \
spf_fo_b03_t02_np33 \
spf_fo_b03_t03_np33 \
spf_fo_b04_t01_np33 \
spf_fo_b04_t02_np33 \
spf_fo_b05_t01_np33 \
spf_fo_b05_t02_np33 \
spf_fo_b05_t03_np33 \
spf_fo_b05_t04_np33 \
> spf_fo_np_rep.csv && \
intensity.r spf_fo_np_rep.csv "SPF vs GF FO Non Productive Repetitions"

somatic surveyPrevalence spf_fo_np33 > spf_fo_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_np33 > gf_fo_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_p33 > spf_fo_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_p33 > gf_fo_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_npn133 > spf_fo_npn133_vj_prevalence.csv
somatic surveyPrevalence spf_fo_npn033 > spf_fo_npn033_vj_prevalence.csv
somatic surveyPrevalence spf_fo_pn133 > spf_fo_pn133_vj_prevalence.csv
somatic surveyPrevalence spf_fo_pn033 > spf_fo_pn033_vj_prevalence.csv

somatic surveyPrevalence spf_fo_b01_t01_np33 > spf_fo_b01_t01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b01_t01_p33 > spf_fo_b01_t01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b01_t02_np33 > spf_fo_b01_t02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b01_t02_p33 > spf_fo_b01_t02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b01_t03_np33 > spf_fo_b01_t03_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b01_t03_p33 > spf_fo_b01_t03_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b02_t01_np33 > spf_fo_b02_t01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b02_t01_p33 > spf_fo_b02_t01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b02_t02_np33 > spf_fo_b02_t02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b02_t02_p33 > spf_fo_b02_t02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t01_np33 > spf_fo_b03_t01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t01_p33 > spf_fo_b03_t01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t02_np33 > spf_fo_b03_t02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t02_p33 > spf_fo_b03_t02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t03_np33 > spf_fo_b03_t03_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b03_t03_p33 > spf_fo_b03_t03_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b04_t01_np33 > spf_fo_b04_t01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b04_t01_p33 > spf_fo_b04_t01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b04_t02_np33 > spf_fo_b04_t02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b04_t02_p33 > spf_fo_b04_t02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t01_np33 > spf_fo_b05_t01_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t01_p33 > spf_fo_b05_t01_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t02_np33 > spf_fo_b05_t02_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t02_p33 > spf_fo_b05_t02_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t03_np33 > spf_fo_b05_t03_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t03_p33 > spf_fo_b05_t03_p33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t04_np33 > spf_fo_b05_t04_np33_vj_prevalence.csv
somatic surveyPrevalence spf_fo_b05_t04_p33 > spf_fo_b05_t04_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t01_np33 > gf_fo_b02_t01_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t01_p33 > gf_fo_b02_t01_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t02_np33 > gf_fo_b02_t02_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t02_p33 > gf_fo_b02_t02_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t03_np33 > gf_fo_b02_t03_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b02_t03_p33 > gf_fo_b02_t03_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t01_np33 > gf_fo_b03_t01_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t01_p33 > gf_fo_b03_t01_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t02_np33 > gf_fo_b03_t02_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t02_p33 > gf_fo_b03_t02_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t03_np33 > gf_fo_b03_t03_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t03_p33 > gf_fo_b03_t03_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t04_np33 > gf_fo_b03_t04_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b03_t04_p33 > gf_fo_b03_t04_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t01_np33 > gf_fo_b04_t01_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t01_p33 > gf_fo_b04_t01_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t02_np33 > gf_fo_b04_t02_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t02_p33 > gf_fo_b04_t02_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t03_np33 > gf_fo_b04_t03_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t03_p33 > gf_fo_b04_t03_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t04_np33 > gf_fo_b04_t04_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b04_t04_p33 > gf_fo_b04_t04_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t01_np33 > gf_fo_b05_t01_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t01_p33 > gf_fo_b05_t01_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t02_np33 > gf_fo_b05_t02_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t02_p33 > gf_fo_b05_t02_p33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t03_np33 > gf_fo_b05_t03_np33_vj_prevalence.csv
somatic surveyPrevalence gf_fo_b05_t03_p33 > gf_fo_b05_t03_p33_vj_prevalence.csv
