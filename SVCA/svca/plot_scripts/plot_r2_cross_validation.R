library(ggplot2)
library(reshape2)
library(gplots)
library(plyr)
library(pheatmap)

working_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/results/IMC_cv////res/'
plot_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/results/IMC_cv//plots/r2/'

annotations = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/patient_data/clinical_data.csv'
# annotations=NA
colours = c('#004D7F'  , '#017100', '#FF9300' )

db_file = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/interactions_db/protein.csv'
db = read.csv(db_file)

col_env = colorRampPalette(c('white', '#FF9300'))(n=100)

w=1500
h=800
base_size = 15
image_dir = working_dir
no_int = FALSE  # NULL, FALSE
averaged=TRUE

select_prots = 1:20


read_pred = function(image_dir){
  
  # get image name
  tmp= strsplit(image_dir, split = '/')[[1]]
  image_name = tmp[length(tmp)]
  print(image_name)
  
  if('results' %in% list.files(image_dir)) image_dir = paste(image_dir, 'results', sep='/')
  
  # get protein names
  tmp = list.files(image_dir, full.names = FALSE)
  tmp = tmp[grep('pred.txt', tmp)]
  protein_names = unlist(lapply(strsplit(tmp, '_'), function(x) return(x[1])))
  
  # get effect name
  tmp = list.files(image_dir, full.names = FALSE)
  tmp = tmp[grep('pred.txt', tmp)]
  effect_names = unlist(lapply(strsplit(tmp, '_'), function(x) return(x[3])))
  
  # get file names
  tmp = list.files(image_dir, full.names = TRUE)
  file_names = tmp[grep('pred.txt', tmp)]
  
  # read results
  tmp = lapply(file_names, function(x) read.table(x, header=TRUE))
  for(i in 1:length(tmp)){
    tmp[[i]] = cbind(effect_names[i], tmp[[i]])
    tmp[[i]] = cbind(protein_names[i], tmp[[i]])
    colnames(tmp[[i]]) = c('protein', 'effect', 'Truth', 'Prediction')
  }
  all_res = do.call(rbind, tmp)
  
  all_res = cbind(image_name, all_res)
  colnames(all_res)[1] = 'image'
  
  return(all_res)
}

calc_r2 = function(truth, pred){
  SS = sum((truth - pred)**2)
  VAR = sum((truth - mean(truth))**2)
  
  return (1 - SS/VAR)
}

calc_all_r2 = function(all_preds_cat){
  all_r2 = ddply(all_preds_cat,  c('image', 'protein', 'effect'), summarise, r2 = calc_r2(Truth, Prediction))
}

scaleFUN <- function(x) sprintf("%.2f", x)

plot_sum_r2 = function(all_r2){
  all_r2[is.infinite(all_r2[,'r2']), 'r2']=0
  all_r2[all_r2[,'r2']<0, 'r2']=0
  
  # select top 20 proteins 
  # all_r2_filt = all_r2[all_r2[,'protein'] %in% res_sum_all_tmp[,'protein'],]
  all_r2_filt = all_r2
  
  sum_r2 = ddply(all_r2_filt, c('effect'), summarise, mean_r2 = mean(r2), sd_r2= sd(r2))
  
  tmp = as.vector(sum_r2[,'effect'])
  tmp[tmp == "local"] = '+Environmental effect'
  tmp[tmp == "env"] = '+Cell-cell interactions'
  tmp[tmp == "intrinsic"] = 'Intrinsic effect'
  sum_r2$effect = as.factor(tmp)
  
  sum_r2[,'effect'] = factor(sum_r2[,'effect'], levels= c('Intrinsic effect', '+Environmental effect', '+Cell-cell interactions'))
  
  sd_error_bars = aes(ymax = mean_r2 + sd_r2, ymin=mean_r2 - sd_r2)
  
  p = ggplot(sum_r2, aes(x=effect, y = mean_r2, fill=effect))+ #, alpha=annotation
    geom_bar(stat="identity", colour='black', position = position_dodge(), alpha=0.75)+
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 15, angle = 0, hjust = 0, colour = "black"),legend.title = element_blank())+
    labs(x = "", y = "R2")+
    scale_fill_manual( values = colours)+ #labels = c("intrinsic", '+environmental', '+cell-cell \n interactions'),
    # scale_fill_manual(values = colours)+
    theme(legend.position="bottom")+
    theme(legend.text=element_text(size=base_size))+
    geom_errorbar(sd_error_bars, width=.3,position=position_dodge(width=.9))+
    scale_y_continuous(labels=scaleFUN)+
    guides(fill=guide_legend(nrow=3,byrow=TRUE))
  # theme(legend.margin=margin(t=0, r=0, b=-0.5, l=0, unit="cm"))
  # guides(fill=guide_legend(nrow=2,byrow=TRUE))
  # theme(legend.position = 'none')
  p  
  ggsave(paste(plot_dir, '/r2_average_mean.pdf', sep = '/'),p, device = 'pdf', width = 3.3, height = 5)
}

all_images = list.files(working_dir, full.names = TRUE)
all_preds = lapply(all_images, read_pred)

all_preds_cat = do.call('rbind', all_preds)
all_r2 = calc_all_r2(all_preds_cat)






