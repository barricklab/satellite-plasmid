require(ggplot2)
require(tidyr)
require(dplyr)

#Run repeats.sh firs to generate MUMmer output for this script to use

# Load data from the three different MUMmer output files
gfp_del_1 = read.table("gfp-del-1.coords.tab", skip=3, sep="\t", row.names=NULL)
colnames(gfp_del_1) = c("S1", "E2", "S2", "S3", "LEN1", "LEN2", "IDY", "SEQ1", "SEQ2")
gfp_del_1$type="deletion"

gfp_del_2 = read.table("gfp-del-2.coords.tab", skip=3, sep="\t", row.names=NULL)
colnames(gfp_del_2) = c("S1", "E2", "S2", "S3", "LEN1", "LEN2", "IDY", "SEQ1", "SEQ2")
gfp_del_2$type="deletion"

gfp_del_3 = read.table("gfp-del-3.coords.tab", skip=3, sep="\t", row.names=NULL)
colnames(gfp_del_3) = c("S1", "E2", "S2", "S3", "LEN1", "LEN2", "IDY", "SEQ1", "SEQ2")
gfp_del_3$type="deletion"

satellite = read.table("satellite.coords.tab", skip=3, sep="\t", row.names=NULL)
colnames(satellite) = c("S1", "E2", "S2", "S3", "LEN1", "LEN2", "IDY", "SEQ1", "SEQ2")
satellite$type="satellite"

full_data = gfp_del_1 %>% rbind(gfp_del_2) %>% rbind(gfp_del_3) %>% rbind(satellite)

ggplot(full_data, aes(LEN1,..density.., color=type)) + geom_freqpoly(binwidth=1)
ggsave("repeats.plot.pdf")

h = hist((full_data %>% filter(type=="deletion"))$LEN1, breaks=c(6:20)+0.5, plot=F)
output_data_1 = data.frame(length = (h$breaks + 0.5)[1:length(h$counts)], count = h$counts, density = h$density, type="deletion")

h = hist((full_data %>% filter(type=="satellite"))$LEN1, breaks=c(6:20)+0.5, plot=F)
output_data_2 = data.frame(length = (h$breaks + 0.5)[1:length(h$counts)], count = h$counts, density = h$density, type="satellite")

output_data = output_data_1 %>% rbind(output_data_2)
write.csv(output_data, "repeats.counts.csv")

