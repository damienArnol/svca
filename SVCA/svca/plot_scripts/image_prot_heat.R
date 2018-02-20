library(ggplot2)
library(reshape2)
library(gplots)
library(plyr)
library(pheatmap)

working_dir = '/tmp/test_svca/'
plot_dir = '/tmp/plots/'

# colours = c('#929292', '#004D7F'  , '#017100','#B51700', '#FF9300' )
colours = c('#929292', '#004D7F'  , '#017100', '#FF9300' )

col_env = colorRampPalette(c('white', '#FF9300'))(n=100)

w=1500
h=800
base_size = 15
image_dir = working_dir
noEnv = FALSE

read_effects = function(image_dir){
  
  if('results' %in% list.files(image_dir)) image_dir = paste(image_dir, 'results', sep='/')
  
  # get protein names
  tmp = list.files(image_dir, full.names = FALSE)
  
  # TODO removes
  tmp = tmp[grep('effect', tmp)] 
  tmp = tmp[grep('interactions', tmp)]

  
  ########
  protein_names = unlist(lapply(strsplit(tmp, '_'), function(x) return(x[1])))
  # print(protein_names)
  
  # get file names
  file_names = list.files(image_dir, full.names = TRUE)
  file_names = file_names[grep('effect', file_names)]
  file_names = file_names[grep('interactions', file_names)] 
  ########
  
  # read results
  tmp = lapply(file_names, function(x) read.table(x, header=TRUE))
  all_res = do.call(rbind, tmp)
  all_res = cbind(protein_names, all_res)
  ix = 2:dim(all_res)[2]
  all_res[,ix] = all_res[,ix]/rowSums(all_res[,ix])
  colnames(all_res)[1] = 'protein'
  
  # summarise results
  all_res_df = melt(all_res, id.vars = 'protein')
  colnames(all_res_df) = c('protein', 'effect', 'value')
  res_sum = ddply(all_res_df, c('protein', 'effect'), summarise, mean_value=mean(value, na.rm=TRUE))
  res_sum$effect <- factor(res_sum$effect, levels=c('noise',"intrinsic", "environmental", 'interactions'))
  
  return(res_sum)
  
}


plot_sig = function(working_dir, plot_dir){
  all_images = list.files(working_dir, full.names = TRUE)
  all_images_short = list.files(working_dir, full.names = FALSE)
  
  all_res = lapply(all_images, read_effects)
  for(i in 1:length(all_res)){
    colnames(all_res[[i]])[3] = all_images_short[i]
  }
  
  all_res_merged = all_res[[1]]
  if (length(all_res) > 1){
    for(i in 2:length(all_res)){# length(all_res)
      all_res_merged = merge(all_res_merged, all_res[[i]], by=c('protein', 'effect'), all=TRUE)  
    }
  }
  
  ###############################################################################################
  # renaming
  ###############################################################################################
  all_res_merged = melt(all_res_merged)
  
  tmp = as.vector(all_res_merged$effect)
  tmp[tmp =='environmental'] = 'Environmental effect'
  tmp[tmp =='interactions'] = 'Cell-cell interactions'
  tmp[tmp =='noise'] = 'Residual noise'
  tmp[tmp =='intrinsic'] = 'Intrinsic effect'
  all_res_merged$effect = tmp
  
  all_res_merged$effect <- factor(all_res_merged$effect, levels=c("Residual noise", "Intrinsic effect", 'Environmental effect', 'Cell-cell interactions'))
  
  ###############################################################################################
  # do violin plot
  p=ggplot(all_res_merged, aes(x=effect, y=value, fill=effect), alpha=0.75)+
    geom_violin() +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(size = 12, angle = 330, hjust = 0, colour = "black"),axis.text.y = element_blank(),legend.title = element_blank())+
    labs(x = "", y = "")+
    scale_fill_manual(values = colours)+
    theme(legend.position="none")+
    coord_flip()
  p
  ggsave(paste(plot_dir, '/violin_plots.pdf', sep = '/'),p, device = 'pdf', width = 2, height = 4)
  
  
  ###############################################################################################
  p = ggplot(all_res_merged, aes(x = protein, y = value, fill=effect))+
    geom_boxplot()+
    scale_alpha_discrete(range=c(0.3, 0.8))+
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "black"),axis.text.y = element_text(size = base_size, angle = 0, hjust = 0, colour = "black"))+
    labs(x = "", y = "variance explained")+
    scale_fill_manual(values = colours)+
    theme(legend.position="bottom")+
    theme(legend.text=element_text(size=base_size))
    # guides(fill=guide_legend(nrow=2,byrow=TRUE))
  # theme(legend.position = 'none')
  ggsave(paste(plot_dir, '/boxplots.png', sep = '/'),p, device = 'png', width = 15, height = 5)

    res_sum = all_res_merged[,c(1, 2, 4)]
    res_sum = ddply(res_sum, c('protein', 'effect'), summarise, variance=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE))
    res_un2 = dcast(res_sum[,-4], protein~effect)
    res_un2$cc= res_un2[,'Cell-cell interactions']
    tmp  = melt(res_un2, id.vars = c('protein', 'cc'))
    colnames(tmp)= c('protein', 'cc', 'effect', 'value')
    res_sum_all = merge(tmp, res_sum, by = c('protein', 'effect'))

  
  # select top 20
  sel = sort(unique(res_sum_all$env), decreasing = T)[1:20]
  # res_sum_all_tmp = res_sum_all[res_sum_all$env %in% sel,]
  res_sum_all_tmp = res_sum_all
  p = ggplot(res_sum_all_tmp, aes(x = reorder(protein, -cc), y = variance, fill=effect))+
    geom_bar(stat="identity", colour='black', alpha=0.75)+
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(size = 15, angle = 330, hjust = 0, colour = "black"),axis.text.y = element_text(size = 15, angle = 0, hjust = 0, colour = "black"),legend.title = element_blank())+
    labs(x = "", y = "Variance explained")+
    scale_fill_manual(values = colours)+
    theme(legend.position="bottom")+
    theme(legend.text=element_text(size=20))
  p
  
  ggsave(paste(plot_dir, '/bars.pdf', sep = '/'), p, device = 'pdf', width = 13, height = 5)

  
}

plot_sig(working_dir, plot_dir)

