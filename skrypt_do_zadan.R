#### Zadanie 1:
library(stringr)

co_druga_wielka <- function(input_string) {
  if (nchar(input_string)>0) { # Checking whether input is not empty
    input_string <- tolower(input_string) # Lower-casing all letters in the string
    str1 <- str2 <- input_string # Create two output variables
    odd <- seq(1, nchar(input_string), 2) # Create indices for the odd letters
    even <- seq(2, nchar(input_string), 2) # Create indices for the even letters
    str1 <- paste(str_sub(str1, odd, odd), collapse='') # Selecting and concatenating letters
    str2 <- paste(str_sub(str2, even, even), collapse='')
    result_list <- list(str1, str2) # Merging into a result list
    result_list # Return the list
  }
  else { print('Empty input') }
}

x <- c('abcdef')
co_druga_wielka(x)
x <- c('abcDEF')
co_druga_wielka(x)
x <- c('')
co_druga_wielka(x)

#### Zadanie 2
liczba_unikalnych_liter <- function(string_vector) {
  if (nchar(string_vector)>0) {
    string_vector <- tolower(string_vector) # Zamiana liter na male
    chars <- str_split(string_vector, '') # Podzielenie stringu na pojedyncze litery
    char_count <- data.frame(table(chars)) # Zliczenie poszczegolnych liter
    char_count <- char_count[char_count$Freq>=2,] # Wybranie tych co maja conajmniej 2 wystapienia
    n_chars <- nrow(char_count) # Zliczenie ile takich jest
    n_chars # Zwrocenie wyniku
  }
  else {
    n_chars <- 0
    n_chars
  }
}

liczba_unikalnych_liter('ABBA')
liczba_unikalnych_liter('aBcbA')
liczba_unikalnych_liter('RabarbArka')
liczba_unikalnych_liter('')
liczba_unikalnych_liter('a')

#### Zadanie 3  
# Subsetting the data for faster manipulations - running a command in the terminal:
# head -n -1000 CPCT02220079.annotated.processed.vcf > test.vcf (from R it is also possible with system() but for some reason it doesn't work...)

library(data.table)
# Reading-in the test data
test_data <- read.csv('test.vcf', skip = 402, sep='\t', stringsAsFactors = F)
full_data <- fread('CPCT02220079.annotated.processed.vcf', skip = 402, sep='\t', stringsAsFactors = F) # fread() is optimised for speed
colnames(test_data)[1] <- 'CHROM'
colnames(full_data)[1] <- 'CHROM'
data <- full_data
subset_chr <- 12 # Subsetting the data:
subset_from <- 112204691
subset_to <- 112247789
data_subset <- data[data$CHROM==subset_chr,]# selecting chromosome
data_subset <- data_subset[data_subset$POS>=subset_from & data_subset$POS<= subset_to,] # selecting position
write.csv(data_subset, paste0('vcf_subset_', subset_chr, '_', subset_from, '_', subset_to, '.vcf'))

#### Zadanie 4 
## Zakladam, ze na wykresach maja byc dlugosci insercji i delecji razem, a nie na osobnych wykresach.
library(ggplot2)
library(dplyr)
library(ggpubr)

data <- full_data
# Counting the insertions and eletions lengths:
data$del_len <- data$ins_len <- NA 
data$del_len <- nchar(data$REF)-1  
data$ins_len <- nchar(data$ALT)-1 
data$indel_len <- data$ins_len + data$del_len 
chromosomes <- unique(data$CHROM)

# Making graphs:
# chr<-'8'
draw_histograms <- function(chr) {
  data_for_plot <- data %>% 
    filter(indel_len>0) %>% 
    filter(CHROM==chr)
  ggplot(data_for_plot, aes(indel_len)) +
    geom_histogram(binwidth = 10) +
    xlab('Indel length') +
    ylab('Indels count log10') + 
    ggtitle(paste0('Chromosome ', chr)) +
    scale_y_log10(labels=function(x) format(x, big.mark = " ", scientific = FALSE))
}
h1 <- draw_histograms(1)
h2 <- draw_histograms(2)
h3 <- draw_histograms(3)
h4 <- draw_histograms(4)
h5 <- draw_histograms(5)
h6 <- draw_histograms(6)
h7 <- draw_histograms(7)
h8 <- draw_histograms(8)
h9 <- draw_histograms(9)
h10 <- draw_histograms(10)
h11 <- draw_histograms(11)
h12 <- draw_histograms(12)
h13 <- draw_histograms(13)
h14 <- draw_histograms(14)
h15 <- draw_histograms(15)
h16 <- draw_histograms(16)
h17 <- draw_histograms(17)
h18 <- draw_histograms(18)
h19 <- draw_histograms(19)
h20 <- draw_histograms(20)
h21 <- draw_histograms(21)
h22 <- draw_histograms(22)
hX <- draw_histograms('X')
hY <- draw_histograms('Y')
hMT <- draw_histograms('MT')

# Drawing all together:
gg <- ggarrange(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, hX, hY, hMT,
            nrow = 5, ncol = 5)
ggsave('Indels.png', gg, height = 15, width=15)

#### Zadanie 5:

library(tidyr)
library(stringr)

#data <- test_data
data <- full_data
data <- data %>% 
  filter(FILTER=='PASS') %>% # Filtering for PASS in FILTER
  filter(grepl('AF=0.5', INFO)) %>% # Filtering for heterozygotes
  filter(grepl('GoNLv5_AF', INFO)) # Filtering for GoNLv5_AF present

data$GoNLv5_AF <- str_split(data$INFO, pattern=';GoNLv5_AF=') %>% 
    sapply( "[", 2 ) 
data$GoNLv5_AF <- str_split(data$GoNLv5_AF, pattern=';GoNLv5_AN=') %>% 
  sapply( "[", 1 ) %>% as.numeric()

data <- data %>% filter(GoNLv5_AF<0.01)  # Filtering for GoNLv5_AF <0.01
#write.csv(x = data, 'data_with_GoNLv5_AF_less_than_0_01.csv')
  
#### Zadanie 6:

#data <- test_data
data <- full_data
data <- data %>% # Selecting the DP - depth of coverage values:
  separate(CPCT02220079R, c('GT', 'AD', 'DP'), sep=':',extra='drop') %>%
  select('CHROM', 'DP') %>% 
  filter(!is.na(as.numeric(CHROM))) # selecting only autosomal chromosomes
data[is.na(data$DP)] <- 0 # Replacing '.' with 0
data$DP <- as.numeric(data$DP)

data_for_plot <- data %>% # Counting the means per chromosome
  group_by(CHROM) %>% 
  summarise(mean = mean(DP))

data_for_plot$CHROM <- factor(data_for_plot$CHROM, levels=sort(as.numeric(data_for_plot$CHROM))) # Sorting the x axis.

ggplot(data_for_plot, aes(CHROM, mean)) + # Creating a plot.
  geom_bar(stat='identity') +
  ylab('Mean depth of coverage') +
  xlab('Chromosomes')
ggsave('Sequencing_depth_per_chromosome.png')  
  
    