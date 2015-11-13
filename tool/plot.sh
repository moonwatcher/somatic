#!/bin/zsh

SESSION_NAME="somatic_plot"
START_DIRECTORY="/Users/lg/thesis/plot"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n etc -d
new-window -t"$SESSION_NAME:1" -n one "\
somatic -v debug surveyPlot spf_b1a_np33         gf_b1a_np33
somatic -v debug surveyPlot spf_b1a_np33         spf_iftl_np33
somatic -v debug surveyPlot spf_b1a_npn033       spf_iftl_np33
somatic -v debug surveyPlot spf_b1a_npn133       spf_b1a_npn033
somatic -v debug surveyPlot spf_b1a_npn133       spf_iftl_np33
somatic -v debug surveyPlot spf_b1a_p33          gf_b1a_p33
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:2" -n two "\
somatic -v debug surveyPlot spf_b1a_p33          spf_iftl_p33
somatic -v debug surveyPlot spf_b1a_pn033        spf_iftl_p33
somatic -v debug surveyPlot spf_b1a_pn133        spf_b1a_pn033
somatic -v debug surveyPlot spf_b1a_pn133        spf_iftl_p33
somatic -v debug surveyPlot spf_fo_np33          gf_fo_np33
somatic -v debug surveyPlot spf_fo_np33          spf_b1a_np33
somatic -v debug surveyPlot spf_fo_np33          spf_b1a_npn033
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:3" -n three "\
somatic -v debug surveyPlot spf_fo_np33          spf_iftl_np33
somatic -v debug surveyPlot spf_fo_p33           gf_fo_p33
somatic -v debug surveyPlot spf_fo_p33           spf_b1a_p33
somatic -v debug surveyPlot spf_fo_p33           spf_b1a_pn033
somatic -v debug surveyPlot spf_fo_p33           spf_iftl_p33
somatic -v debug surveyPlot spf_iabm_np33
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:4" -n four "\
somatic -v debug surveyPlot spf_iabm_p33
somatic -v debug surveyPlot spf_iaspl_np33
somatic -v debug surveyPlot spf_iaspl_p33
somatic -v debug surveyPlot spf_iftl_np33
somatic -v debug surveyPlot spf_iftl_p33
somatic -v debug surveyPlot spf_mz_np33          gfmz_np33
somatic -v debug surveyPlot spf_mz_p33           gfmz_p33
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:5" -n five "\
somatic -v debug surveyPlot spf_np33             gf_np33
somatic -v debug surveyPlot spf_p33              gf_p33
somatic -v debug surveyPlot spf_preb_np33        spf_b1a_np33
somatic -v debug surveyPlot spf_preb_np33        spf_fo_np33
somatic -v debug surveyPlot spf_preb_np33        spf_mz_np33
somatic -v debug surveyPlot spf_preb_p33         spf_b1a_p33
somatic -v debug surveyPlot spf_preb_p33         spf_fo_p33
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:6" -n six "\
somatic -v debug surveyPlot spf_preb_p33         spf_mz_p33
somatic -v debug surveyPlot spf_prebabm_np33
somatic -v debug surveyPlot spf_prebabm_p33
somatic -v debug surveyPlot spf_prebftl_np33
somatic -v debug surveyPlot spf_prebftl_p33
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"


