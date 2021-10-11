# Clustering overview

##### A) LOAD PACKAGES AND SET PARAMETERS ########
library(ggplot2)
library(optparse)

# (Semi-) hardcoded parameters for plotting
params <- vector(mode = "list")
params$plotwidth = 7
params$plotheight = 5
params$plotunit = 'in'
params$debug = F

##### B) DIGEST USER INPUT WITH OPTPARSE#######
if (!params$debug){
  option_list <- list(
    make_option(c("-i", "--input_file"), type = "character", default = NULL,
                help = "", metavar = "character"),
    make_option(c("-o", "--output_plot"), type = "character", default = NULL,
                help = "", metavar = "character")
  )


opt <- parse_args(OptionParser(option_list = option_list))

infile_link    = opt$input_file 
outplot_link  = opt$output_plot
} else if (params$debug){
  # For easy running in debug mode
  infile_link = "../data/n_fields.txt"
  output_plot = "../plots/fooplot1.pdf"
}


clusters = read.table(infile_link, sep=' ')

p = ggplot(clusters) + geom_point(aes(x=1:length(V1), y=sort(V1, decreasing=T))) + 
  geom_line(aes(x=1:length(V1), y=sort(V1, decreasing=T))) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x='# Cluster', y='Reads belonging to cluster', 
       title= paste0('CARNAC-LR: Reads and clusters \nSum or reads: ',
                     sum(clusters),'\nSingletons: ', sum(clusters$V1==1),
                     ' clusters. >= 2 reads: ', sum(clusters$V1>1), 
                     ' clusters'))

ggsave(p, file=outplot_link, device='pdf', width = params$plotwidth, 
       height=params$plotheight, units=params$plotunit)
