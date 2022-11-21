options(scipen = 999)
library(tidyverse)
library(ggplot2)
test_length <- read.delim("testLength.results",header = TRUE,sep = "\t")
# First figure
pdf(file="bpLength.pdf", width=8, height=4) 
test_length %>% 
  pivot_longer(!Length, names_to = "Reads", values_to = "Counts") %>% 
  mutate(Reads = as.factor(Reads)) %>% 
  mutate(Reads = fct_relevel(.f = Reads, "Filtered", "Uniques", "Otus")) %>% 
	ggplot(aes(x=Length, y=Counts, color=Reads)) +
	geom_line(lwd = 1.2, linetype = 1) +
	theme(legend.position = "bottom") +
	labs(title="Filtered, uniques, OTUs number per bp length", x="Length (bp)",y="Number of reads/OTUs") +
	guides(color=guide_legend(title=NULL))
dev.off()
# Second figure
pdf(file="bpLength-log2.pdf", width=8, height=4) 
test_length %>% 
  pivot_longer(!Length, names_to = "Reads", values_to = "Counts") %>% 
  mutate(Reads = as.factor(Reads)) %>% 
  mutate(Reads = fct_relevel(.f = Reads, "Filtered", "Uniques", "Otus")) %>% 
	ggplot(aes(x=Length, y=Counts, color=Reads)) +
	geom_line(lwd = 1.2, linetype = 1) +
	theme(legend.position = "bottom") +
	scale_y_continuous(trans="log2") +
	labs(title="Filtered, uniques, OTUs number per bp length", x="Length (bp)",y="Log2 number of reads/OTUs") +
	guides(color=guide_legend(title=NULL))
dev.off()
pdf(file="bpLength-faceted.pdf", width=12, height=4) 
test_length %>% 
  pivot_longer(!Length, names_to = "Reads", values_to = "Counts") %>% 
  mutate(Reads = as.factor(Reads)) %>% 
  mutate(Reads = fct_relevel(.f = Reads, "Filtered", "Uniques", "Otus")) %>% 
	ggplot(aes(x=Length, y=Counts, color=Reads)) +
	geom_line(lwd = 1.2, linetype = 1) +
	theme(legend.position = "bottom") +
	facet_wrap(~Reads, scales = "free") +
	labs(title="Filtered, uniques, OTUs number per bp length", x="Length (bp)",y="number of reads/OTUs") +
	guides(color=guide_legend(title=NULL))
dev.off()
