library(readxl)
library(stringr)
library(seqinr)
library(purrr)


freq_HC_mutations = read.table('<Convergent heavy chain mutation file>', header = FALSE)
freq_HC_mut_df = data.frame(Mutations = freq_HC_mutations$V1, Chain = rep('HC',length(freq_HC_mutations$V1)))
freq_KC_mutations = read.table('<Convergent light chain mutation file>', header = FALSE)
freq_KC_mut_df = data.frame(Mutations = freq_KC_mutations$V1, Chain = rep('KC',length(freq_KC_mutations$V1)))

freq_mut_df = rbind(freq_HC_mut_df, freq_KC_mut_df)


lineage_fasta_file = read.fasta(file = '<Heavy or light chain fasta file with M15-like sequences from donor>', seqtype = "AA")
lineage_headers = names(lineage_fasta_file)
lineage_headers = lineage_headers[lineage_headers != "UCA"]

lineage_hc_file = read.table('<Heavy chain mutation file from the donor>', header = F)
lineage_hc_df = cbind(rep('<Donor ID>',length(lineage_hc_file$V1)), lineage_hc_file, rep('HC', length(lineage_hc_file$V1)))
colnames(lineage_hc_df) = c('Subject', 'Mutation', 'Sequences', 'Chain')

lineage_kc_file = read.table('<Light chain mutation file from the donor>', header = F)
lineage_kc_df = cbind(rep('Donor-1',length(lineage_kc_file$V1)), lineage_kc_file, rep('KC', length(lineage_kc_file$V1)))
colnames(lineage_kc_df) = c('Subject', 'Mutation', 'Sequences', 'Chain')

lineage_mut_df = rbind(lineage_hc_df, lineage_kc_df)

lineage_mut_linkage_OR_matrix = matrix(data = NA, nrow = length(freq_mut_df$Mutations), ncol = length(freq_mut_df$Mutations))


count = 1

lineage_mutation_1 = c()
lineage_mutation_2 = c()
lineage_odds_ratio = c()
lineage_mutation_1_chain = c()
lineage_mutation_2_chain = c()


for(i in 1:length(freq_mut_df$Mutations)){
  
  for(j in 1:length(freq_mut_df$Mutations)){
    
    lineage_mutation_1_chain[count] = freq_mut_df$Chain[i]
    lineage_mutation_2_chain[count] = freq_mut_df$Chain[j]
    
    if((freq_mut_df$Mutations[i] %in% lineage_mut_df$Mutation) && (freq_mut_df$Mutations[j] %in% lineage_mut_df$Mutation) && i<j){
      
      
      contingency_matrix = matrix(data = 0, nrow = 2, ncol = 2)
      mut1_df = lineage_mut_df[lineage_mut_df$Mutation == freq_mut_df$Mutations[i],]
      mut1_seq = as_vector(strsplit(mut1_df$Sequences[1],"," ))
      mut2_df = lineage_mut_df[lineage_mut_df$Mutation == freq_mut_df$Mutations[j],]
      mut2_seq = as_vector(strsplit(mut2_df$Sequences[1],"," ))
      not_mut1_seq = lineage_headers[!(lineage_headers %in% mut1_seq)]
      not_mut2_seq = lineage_headers[!(lineage_headers %in% mut2_seq)]
      
      contingency_matrix[1,1] = length(intersect(mut1_seq, mut2_seq))
      contingency_matrix[1,2] = length(intersect(mut1_seq, not_mut2_seq))
      contingency_matrix[2,1] = length(intersect(not_mut1_seq, mut2_seq))
      contingency_matrix[2,2] = length(intersect(not_mut1_seq, not_mut2_seq))
      
      rownames(contingency_matrix) = c(freq_mut_df$Mutations[i], paste("Not_",freq_mut_df$Mutations[i], sep = ''))
      colnames(contingency_matrix) = c(freq_mut_df$Mutations[j], paste("Not_",freq_mut_df$Mutations[j], sep = ''))
      
      f_test = fisher.test(contingency_matrix, alternative = "greater")
      

      lineage_mut_linkage_OR_matrix[i,j] = f_test$estimate[[1]]
      
      lineage_mutation_1[count] = freq_mut_df$Mutations[i]
      lineage_mutation_2[count] = freq_mut_df$Mutations[j]
      
      lineage_odds_ratio[count] = f_test$estimate[[1]]



      
      
     
      count = count + 1
        
    } else{
      
      lineage_mutation_1[count] = freq_mut_df$Mutations[i]
      lineage_mutation_2[count] = freq_mut_df$Mutations[j]
      lineage_odds_ratio[count] = 1
      
      count = count + 1
    }
    
    
    
  }
  
  
  
  
}

rownames(lineage_mut_linkage_OR_matrix) = freq_mut_df$Mutations
colnames(lineage_mut_linkage_OR_matrix) = freq_mut_df$Mutations


write.csv(lineage_mut_linkage_OR_matrix, '<CSV file with Odds ratios>')

