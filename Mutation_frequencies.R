library(tidyverse)
library(stringr)
library(ggpubr)
library(readxl)


#The following code takes six files as input: M15-like sequences mutation file, Influenza mutation file, SARS-CoV-2 (non-S2) mutation file, Malaria mutation file, sequence count file and convergent SHMs file. The code presented here is for one chain (heavy/light)

#COVID-19

M15_file = read.table('<Path to M15-like sequences mutation file>', header = F)
M15_df = cbind(rep('M15-like_sequences',length(M15_file$V1)), M15_file)
colnames(M15_df) = c('Subject', 'Mutation', 'Sequences')

all_clonotypes_df = M15_df

#Influenza

Influenza_file = read.table('<Path to Influenza HV4-59 or KV3-20 mutation file>', header = F)
Influenza_df = cbind(rep('Influenza',length(Influenza_file$V1)), Influenza_file)
colnames(Influenza_df) = c('Subject', 'Mutation', 'Sequences')



all_clonotypes_Influenza_df = rbind(all_clonotypes_df, Influenza_df)

#CoV

CoV_hc_file = read.table('<Path to SARS-CoV-2 non-S2 HV4-59 or KV3-20 mutation file>', header = F)
CoV_df = cbind(rep('SARS-CoV-2(non-S2)',length(CoV_file$V1)), CoV_file)
colnames(CoV_df) = c('Subject', 'Mutation', 'Sequences')



all_clonotypes_CoV_df = rbind(all_clonotypes_Influenza_df, CoV_df)


#Malaria

Malaria_file = read.table('<Path to Malaria HV4-59 or KV3-20 mutation file>', header = F)
Malaria_df = cbind(rep('Malaria',length(Malaria_file$V1)), Malaria_file)
colnames(Malaria_df) = c('Subject', 'Mutation', 'Sequences')



all_clonotypes_Malaria_df = rbind(all_clonotypes_CoV_df, Malaria_df)



sequence_count_file = read.table('<Text file with the number of sequences in each antigen group>', header = T)

mutations = read.table('<Text file with the convergent mutations>', header = F)


mutation_frequency_matrix = matrix(data = 0, nrow = length(sequence_count_file$Subject), ncol = length(mutations$V1))



count = 1

for(i in 1:length(sequence_count_file$Subject)){
  subject_df = all_clonotypes_Malaria_df[all_clonotypes_Malaria_df$Subject == sequence_count_file$Subject[i],]
  
  
  for(j in 1:length(mutations$V1)){
    
    
    sub_df = subject_df[subject_df$Mutation == mutations$V1[j],]
    
    if(length(sub_df$Subject) == 0){
      
      mutation_frequency_matrix[i,j] = 0
      
      
    }
    
    else{
      
      mutation_frequency_matrix[i,j] = length(as_vector(strsplit(sub_df$Sequences[1], ",")))*100/sequence_count_file$Count[i]
            
    }
    
    count = count + 1
    
    
  }
  
  
}
rownames(mutation_frequency_matrix) = sequence_count_file$Subject
colnames(mutation_frequency_matrix) = hc_mutations$V1

write.csv(mutation_frequency_matrix, '<CSV file convergent mutation frequency matrix>')


