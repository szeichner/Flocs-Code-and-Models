# Analyze statistics for flocs data

data <- read.csv("SettlingVelocities_forR.csv")

# Use t.test to quantify different assertions of difference in the text:

# (1) Flocculated clays+polymers were faster than controls
# control: 1, 12
# flocculated: 2:11, 13:29 
control = data[c(1,12),6]
mean(control)
flocs = data[c(2:11, 13:29),6]
mean(flocs)
t.test(flocs, control, alternative="greater")

# (2) Guar gum was a better flocculant than xanthan
# guar =  2:6, 13:17, 23:26 
# xanthan = 7:11, 18:22,  27:29

guar = data[c(2:6,13:17,  23:26),6]
mean(guar)
xanthan = data[c(7:11, 18:22, 27:29),6]
mean(xanthan)
t.test(guar, xanthan, alternative="greater")

# (3) Guar/smectite combo versus everything else
smectite_guar = data[c(2:6, 23:26),6]
mean(smectite_guar)
everything_else = data[c(7:22, 27:29),6]
mean(everything_else)
t.test(smectite_guar, everything_else, alternative="greater" )
