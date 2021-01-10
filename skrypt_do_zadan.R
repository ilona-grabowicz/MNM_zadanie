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
library(data.table)

# Subsetting the data for faster manipulations - running a command in the terminal:
# head -n -1000 CPCT02220079.annotated.processed.vcf > test.vcf (from R it is also possible with system() but for some reason it doesn't work...)

test_data <- read.csv('test.vcf', skip = 402, sep='\t', stringsAsFactors = F)
  
library(vcfR)
vcf <- read.vcfR('test.vcf', verbose = FALSE )

# read two times the vcf file, first for the columns names, second for the data
tmp_vcf<-readLines("CPCTtest.annotated.processed.vcf")
tmp_vcf_data<-read.table("CPCTtest.annotated.processed.vcf", stringsAsFactors = FALSE)

# filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(tmp_vcf_data)<-vcf_names
vcf <- read.csv('CPCTtest.annotated.processed.vcf', skip=401, sep='\t', stringsAsFactors = F)
  
  
  
  
  
  
  



