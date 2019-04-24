# plots
xterm -e "sleep 0;   python exsan.py -b batch_viii_all_decay.txt" &
xterm -e "sleep 60;  python exsan.py -b batch_viii_all_mg.txt" &
xterm -e "sleep 120; python exsan.py -b batch_viii_all_pointwise.txt" &
# analysis
xterm -e "sleep 180; python exsan.py -b batch_viii_analysis_decay.txt" &
xterm -e "sleep 240; python exsan.py -b batch_viii_analysis_mg.txt" &
xterm -e "sleep 300; python exsan.py -b batch_viii_analysis_pointwise.txt" &

# combinations
# xterm -e "sleep 360; python exsan.py -b batch_viii_multi.txt" &
# xterm -e "sleep 420; python exsan.py -b batch_viii_various.txt" &
