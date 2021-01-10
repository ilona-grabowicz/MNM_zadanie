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


  
  
  
  
  
  
  
  
  
  



