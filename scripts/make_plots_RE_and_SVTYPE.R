# Whoeps, 11th Oct 2021
# Plot RE (Task #2), SVTYPE distribution (Task #4) and INS-RE distribution (Task #5)

##### A) LOAD PACKAGES AND SET PARAMETERS ########
library(ggplot2)
library(dplyr)
library(optparse)

# (Semi-) hardcoded parameters for plotting
params <- vector(mode = "list")
params$plotwidth = 5
params$plotheight = 4
params$plotunit = 'in'
params$debug = F

##### B) DIGEST USER INPUT WITH OPTPARSE#######
if (!params$debug){
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "", metavar = "character"),
  make_option(c("-p", "--output_plot1"), type = "character", default = NULL,
              help = "", metavar = "character"),
  make_option(c("-q", "--output_plot2"), type = "character", default = NULL,
              help = "", metavar = "character"),
  make_option(c("-r", "--output_plot3"), type = "character", default = NULL,
              help = "", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

infile_link    = opt$input_file 
outplot1_link  = opt$output_plot1
outplot2_link  = opt$output_plot2
outplot3_link  = opt$output_plot3
} else if (params$debug){
  # For easy running in debug mode
  infile_link = "../../RE_SVTYPE.tsv"
  output_plot1 = "../plots/fooplot1.pdf"
  output_plot2 = "../plots/fooplot2.pdf"
  output_plot3 = "../plots/fooplot3.pdf"
  
}

###### C) RUN CODE AND MAKE PLOTS ###############

# Load files
RElist = read.table(infile_link)
colnames(RElist) = c('RE','SVTYPE')
  
# Count SV types for nicer plotting of plot #2
RElist2 = RElist %>% 
  count(SVTYPE) %>% 
  mutate(perc = n / nrow(RElist))

# First plot: Histogram of RE, all SVs
p1 = ggplot(RElist) + 
  geom_histogram(aes(x=RE, fill=SVTYPE), bins=100) +
  scale_x_log10() + 
  theme_bw() + 
  labs(x='Reads supporting SV', y='#SVs', title = 'Reads supporting SV')

# Second plot: Barplot of SVTYPES
p2 = ggplot(RElist2) + 
  geom_bar(aes(x=reorder(SVTYPE, -perc), y=n, fill=SVTYPE), stat='identity') + 
  theme_bw() + 
  theme(legend.position = 'none') + 
  scale_y_log10() + 
  labs(x='SV type', y='#SVs', title = 'Frequency of SV types')

# Third plot: Histogram of RE, INS only
p3 = ggplot(RElist[RElist$SVTYPE=='INS',]) + 
  geom_histogram(aes(x=RE), bins=100, fill='#00B8E5') +
  scale_x_log10() + 
  theme_bw() + 
  labs(x='Reads supporting Insertions', y='#INS', title = 'Reads supporting INS')


# Save everything
ggsave(p1, file=outplot1_link, device='pdf', width = params$plotwidth, 
       height=params$plotheight, units=params$plotunit)
ggsave(p2, file=outplot2_link, device='pdf', width = params$plotwidth, 
       height=params$plotheight, units=params$plotunit)
ggsave(p3, file=outplot3_link, device='pdf', width = params$plotwidth, 
       height=params$plotheight, units=params$plotunit)

# Inform user
print('Output plots written.')
