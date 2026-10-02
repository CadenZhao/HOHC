# Load packages ----
set.seed(1)
library(tidyverse)
library(ggsankey)
library(RColorBrewer)
# library(comprehenr)

# plot
plot_sankey <- function (sex, hormone_class, hormone_name, hormone_tissue,
                         receptor_class, receptor_name, receptor_tissue,
                         number_HRtissue, tpm_log1p_1, tpm_log1p_2, gtex_hr, df, vec_receptor){
  ##################################################################
  # get 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  if (sex == 'Male'){
    # define male tissue and male tissue color
    sex_tissue <- c('Testis', 'Prostate')
    sex_tissue_oppose <- c('Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri')
    node_order_sex <- c("Peptide","Lipid-derived","Amino acid-derived","Nuclear HR","Membrane HR",
                        "Testis","Prostate","Thyroid","Stomach","Spleen","Small Intestine","Skin",
                        "Salivary Gland","Pituitary","Pancreas","Nerve","Muscle","Lung","Liver",
                        "Kidney","Heart","Esophagus","Colon","Breast","Brain","Blood Vessel",
                        "Blood","Bladder","Adrenal Gland","Adipose Tissue")
    node_order_sex_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[68:69],col_vector[7:29])
    tissue_addrow <- c("Spleen","Esophagus","Prostate","Muscle","Breast" ,"Bladder")
  } else {
    # define female tissue and female tissue color
    sex_tissue <- c('Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri')
    sex_tissue_oppose <- c('Testis', 'Prostate')
    node_order_sex <- c("Peptide","Lipid-derived","Amino acid-derived","Nuclear HR","Membrane HR",
                        "Vagina","Uterus","Ovary","Fallopian Tube","Cervix Uteri","Thyroid","Stomach",
                        "Spleen","Small Intestine","Skin","Salivary Gland","Pituitary","Pancreas",
                        "Nerve","Muscle","Lung","Liver","Kidney","Heart","Esophagus","Colon","Breast",
                        "Brain","Blood Vessel","Blood","Bladder","Adrenal Gland","Adipose Tissue")
    node_order_sex_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[70:74],col_vector[7:29])
    tissue_addrow <- c("Muscle","Uterus","Vagina","Breast","Spleen","Esophagus","Cervix Uteri","Fallopian Tube","Bladder")
  }
  names(node_order_sex_color) <- node_order_sex
  
  ##################################################################
  # 计算HR表达top的组织
  df_ <- df %>% filter(!Hormone1_tissue %in% sex_tissue_oppose) # 只留性别该有的分泌激素的组织
  gtex_hr_ <- filter(gtex_hr, Description %in% vec_receptor, SEX == sex)
  # 计算最高表达每个HR的组织类型（GTEx数据，30种）
  gtex_hr_ <- gtex_hr_ %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
  gtex_hr_ <- gtex_hr_ %>% rename('Gene_symbol'='Description')
  gtex_hr_ <- gtex_hr_ %>% group_by(Gene_symbol) %>% filter(log1p(TPM) > tpm_log1p_1, log1p(TPM) < tpm_log1p_2, TPM >= sort(TPM, TRUE)[number_HRtissue]) # 此阈值的选定见Fig. S2B
  df_merged <- left_join(df_, gtex_hr_, by='Gene_symbol', relationship = "many-to-many")
  df_merged <- df_merged %>% add_row(Hormone1_tissue = tissue_addrow) # 补全该性别的激素组织
  df_merged <- df_merged %>% rename('Hormone class' = 'Hormone_chemical_classes',
                                    'Hormone name' = 'Hormone_name1',
                                    'Hormone organ' = 'Hormone1_tissue',
                                    'HR organ' = 'SMTS',
                                    'HR type' = 'Receptor_subcellular',
                                    'HR gene name' = 'Gene_symbol')

  # df_merged$`HR type` <- df_merged$`HR type` %>% recode(Membrane = 'Membrane HR', Nucleus = 'Nucleus HR')
  
  ###########################
  # filter condition
  if (hormone_class != 'All') {
    df_sankey <- df_merged %>% filter(`Hormone class` == hormone_class)
  } else {
    df_sankey <- df_merged
  }
  if (!is.null(hormone_name)) {
    df_sankey <- df_sankey %>% filter(`Hormone name` %in% hormone_name)
  }
  if (!is.null(hormone_tissue)) {
    df_sankey <- df_sankey %>% filter(`Hormone organ` %in% hormone_tissue)
  }
  if (receptor_class != 'All') {
    df_sankey <- df_sankey %>% filter(`HR type` == receptor_class)
  }
  if (!is.null(receptor_name)) {
    df_sankey <- df_sankey %>% filter(`HR gene name` %in% receptor_name)
  }
  if (!is.null(receptor_tissue)) {
    df_sankey <- df_sankey %>% filter(`HR organ` == receptor_tissue)
  }
  df_sankey <- df_sankey %>% make_long(`Hormone class`, `Hormone organ`, `HR organ`, `HR type`)
  df_sankey$node <- factor(df_sankey$node, levels = node_order_sex) # 按照性别通用组织排序，然后再排性别特异的组织
  df_sankey$next_node <- factor(df_sankey$next_node, levels = node_order_sex) # 按照性别通用组织排序，然后再排性别特异的组织
  
  if (dim(df_sankey)[1] == 0) {
    print('No results found. You may want to reset the filters or try different options.')
  } else {
    ggplot(df_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
      geom_sankey() +
      geom_sankey_label() +
      theme_sankey() +
      theme(legend.position = "none")+
      scale_fill_manual(values = node_order_sex_color) +
      labs(x = NULL)
  }
}

# download data
filter_data <- function (sex, hormone_class, hormone_name, hormone_tissue,
                         receptor_class, receptor_name, receptor_tissue,
                         number_HRtissue, tpm_log1p_1, tpm_log1p_2, gtex_hr, df, vec_receptor){
  ##################################################################
  # get 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  if (sex == 'Male'){
    # define male tissue and male tissue color
    sex_tissue <- c('Testis', 'Prostate')
    sex_tissue_oppose <- c('Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri')
    node_order_sex <- c("Peptide","Lipid-derived","Amino acid-derived","Nuclear HR","Membrane HR",
                        "Testis","Prostate","Thyroid","Stomach","Spleen","Small Intestine","Skin",
                        "Salivary Gland","Pituitary","Pancreas","Nerve","Muscle","Lung","Liver",
                        "Kidney","Heart","Esophagus","Colon","Breast","Brain","Blood Vessel",
                        "Blood","Bladder","Adrenal Gland","Adipose Tissue")
    node_order_sex_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[68:69],col_vector[7:29])
    tissue_addrow <- c("Spleen","Esophagus","Prostate","Muscle","Breast" ,"Bladder")
  } else {
    # define female tissue and female tissue color
    sex_tissue <- c('Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri')
    sex_tissue_oppose <- c('Testis', 'Prostate')
    node_order_sex <- c("Peptide","Lipid-derived","Amino acid-derived","Nuclear HR","Membrane HR",
                        "Vagina","Uterus","Ovary","Fallopian Tube","Cervix Uteri","Thyroid","Stomach",
                        "Spleen","Small Intestine","Skin","Salivary Gland","Pituitary","Pancreas",
                        "Nerve","Muscle","Lung","Liver","Kidney","Heart","Esophagus","Colon","Breast",
                        "Brain","Blood Vessel","Blood","Bladder","Adrenal Gland","Adipose Tissue")
    node_order_sex_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[70:74],col_vector[7:29])
    tissue_addrow <- c("Muscle","Uterus","Vagina","Breast","Spleen","Esophagus","Cervix Uteri","Fallopian Tube","Bladder")
  }
  names(node_order_sex_color) <- node_order_sex
  
  ##################################################################
  # 计算HR表达top的组织
  df_download <- df %>% filter(!Hormone1_tissue %in% sex_tissue_oppose) # 只留性别该有的分泌激素的组织
  gtex_hr_download <- filter(gtex_hr, Description %in% vec_receptor, SEX == sex)
  # 计算最高表达每个HR的组织类型（GTEx数据，30种）
  gtex_hr_download <- gtex_hr_download %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
  gtex_hr_download <- gtex_hr_download %>% rename('Gene_symbol'='Description')
  gtex_hr_download <- gtex_hr_download %>% group_by(Gene_symbol) %>% filter(log1p(TPM) > tpm_log1p_1, log1p(TPM) < tpm_log1p_2, TPM >= sort(TPM, TRUE)[number_HRtissue]) # 此阈值的选定见Fig. S2B
  df_merged <- left_join(df_download, gtex_hr_download, by='Gene_symbol', relationship = "many-to-many")
  df_merged <- df_merged %>% add_row(Hormone1_tissue = tissue_addrow) # 补全该性别的激素组织
  df_merged <- df_merged %>% rename('Hormone class' = 'Hormone_chemical_classes',
                                    'Hormone name' = 'Hormone_name1',
                                    'Hormone organ' = 'Hormone1_tissue',
                                    'HR organ' = 'SMTS',
                                    'HR type' = 'Receptor_subcellular',
                                    'HR gene name' = 'Gene_symbol')
  
  # df_merged$`HR type` <- df_merged$`HR type` %>% recode(Membrane = 'Membrane HR', Nucleus = 'Nucleus HR')
  
  ###########################
  # filter condition
  if (hormone_class != 'All') {
    df_sankey <- df_merged %>% filter(`Hormone class` == hormone_class)
  } else {
    df_sankey <- df_merged
  }
  if (!is.null(hormone_name)) {
    df_sankey <- df_sankey %>% filter(`Hormone name` %in% hormone_name)
  }
  if (!is.null(hormone_tissue)) {
    df_sankey <- df_sankey %>% filter(`Hormone organ` %in% hormone_tissue)
  }
  if (receptor_class != 'All') {
    df_sankey <- df_sankey %>% filter(`HR type` == receptor_class)
  }
  if (!is.null(receptor_name)) {
    df_sankey <- df_sankey %>% filter(`HR gene name` %in% receptor_name)
  }
  if (!is.null(receptor_tissue)) {
    df_sankey <- df_sankey %>% filter(`HR organ` == receptor_tissue)
  }
  df_sankey %>% drop_na(`HR gene name`)
}