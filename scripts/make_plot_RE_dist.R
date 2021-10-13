# Whoeps, 11th Oct 2021
# Plot RE of long INS

##### A) LOAD PACKAGES AND SET PARAMETERS ########
library(ggplot2)
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
    make_option(c("-o", "--output_plot"), type = "character", default = NULL,
                help = "", metavar = "character")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  infile_link    = opt$input_file 
  outplot_link  = opt$output_plot
} else if (params$debug){
  # For easy running in debug mode
  infile_link = "../../HG002.Hifi.GRCh37.table_RE_largeINS.tsv"
  output_plot = "../plots/fooplot1.pdf"
}

###### C) RUN CODE AND MAKE PLOTS ###############

# Load files
RElist = read.table(infile_link)

# Make plot
p = ggplot(RElist) + 
  geom_histogram(aes(x=V1), color='white',  fill='#00B8E5', bins = 100) +
  #scale_x_log10() + 
  theme_bw() + 
  theme(legend.position = 'none') + 
  labs(x='Reads supporting long INS', y='#INS', title = 'Reads supporting long INS') 
p
# Save plot
ggsave(p, file=outplot_link, device='pdf', width = params$plotwidth, 
       height=params$plotheight, units=params$plotunit)

# Inform user
print('Output plot written.')


