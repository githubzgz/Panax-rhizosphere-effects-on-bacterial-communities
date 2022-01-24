#Phylogenetic analysis of rhizosphere-enriched taxa


#1.Phylogenetic sampling theory
library(picante)
#(1)Generating data
Merge_func <- function(x,y){
  df <- merge(x,y,by='OTU',all.x = T,all.y = T)
  return(df)
}
all_differ <- Reduce(Merge_func,list(
  data.frame(OTU = c(gin_two_up,gin_three_up,gin_four_up,
                     qui_two_up,qui_three_up,qui_four_up,
                     not_two_up,not_three_up,not_four_up) %>% .[!duplicated(.)],
             num = 1),
  data.frame(OTU = c(gin_two_down,gin_three_down,gin_four_down,
                     qui_two_down,qui_three_down,qui_four_down,
                     not_two_down,not_three_down,not_four_down) %>% .[!duplicated(.)],
             num = 1),
  data.frame(OTU = c(gin_two_nosig,gin_three_nosig,gin_four_nosig,
                     qui_two_nosig,qui_three_nosig,qui_four_nosig,
                     not_two_nosig,not_three_nosig,not_four_nosig) %>% .[!duplicated(.)],
             num = 1)))
#gin_two_up:The rhizosphere-enriched taxa in the 2-year-old Panax ginseng field;
#gin_two_down: The rhizosphere-depleted taxa in the 2-year-old Panax ginseng field.
rownames(all_differ) <- all_differ$OTU
all_differ$OTU <- NULL
all_differ[is.na(all_differ)]<- 0
names(all_differ) <- c('up','down','nosig')

all_differ_match <- match.phylo.comm(b_tree,all_differ%>%t)#b_tree: phylogenetic tree

#(2)Calculating the real phylogenetic diversity (PD) and expected PD
all_real_pd <- pd(all_differ_match$comm,all_differ_match$phy)
all_expect_pd <- variance.pd(all_differ_match$phy,upper.bound = F) #variance.pd -- picante package

#(3)Calculated the z-score
up_z_score <- (all_differ_pd[1,1] - all_expect_pd[all_differ_pd[1,2],1])/sqrt(all_expect_pd[all_differ_pd[1,2],2])
down_z_score <- (all_differ_pd[2,1] - all_expect_pd[all_differ_pd[2,2],1])/sqrt(all_expect_pd[all_differ_pd[2,2],2])
other_z_score <- (all_differ_pd[3,1] - all_expect_pd[all_differ_pd[3,2],1])/sqrt(all_expect_pd[all_differ_pd[3,2],2])

#2.Resampling-based unweighted UniFrac distances
my_rarefy_mean <- function(gin_two_up,gin_two_down,gin_two_nosig,gin_two_all,
                           gin_three_up,gin_three_down,gin_three_nosig,gin_three_all,
                           gin_four_up,gin_four_down,gin_four_nosig,gin_four_all,
                           qui_two_up,qui_two_down,qui_two_nosig,qui_two_all,
                           qui_three_up,qui_three_down,qui_three_nosig,qui_three_all,
                           qui_four_up,qui_four_down,qui_four_nosig,qui_four_all,
                           not_two_up,not_two_down,not_two_nosig,not_two_all,
                           not_three_up,not_three_down,not_three_nosig,not_three_all,
                           not_four_up,not_four_down,not_four_nosig,not_four_all,
                           tree,rand=100){
  require(phyloseq)
  require(tidyr)
  require(reshape2)
  require(stringr)
  require(vegan)
  min_num <- min(c(length(gin_two_up),length(gin_two_down),length(gin_two_nosig),length(gin_two_all),
                   length(gin_three_up),length(gin_three_down),length(gin_three_nosig),length(gin_three_all),
                   length(gin_four_up),length(gin_four_down),length(gin_four_nosig),length(gin_four_all),
                   length(qui_two_up),length(qui_two_down),length(qui_two_nosig),length(qui_two_all),
                   length(qui_three_up),length(qui_three_down),length(qui_three_nosig),length(qui_three_all),
                   length(qui_four_up),length(qui_four_down),length(qui_four_nosig),length(qui_four_all),
                   length(not_two_up),length(not_two_down),length(not_two_nosig),length(not_two_all),
                   length(not_three_up),length(not_three_down),length(not_three_nosig),length(not_three_all),
                   length(not_four_up),length(not_four_down),length(not_four_nosig),length(not_four_all)))#È¡×îÐ¡Öµ
  
  uuf_result <- array(NA,dim = c(36,36,rand))
  for(i in 1:rand){
    print(paste(i,'in',rand,sep=''))
    otu <- Reduce(Merge_func,list(data.frame(OTU = gin_two_up[sample(length(gin_two_up),min_num)],gin_two_up=1),
                                  data.frame(OTU = gin_two_down[sample(length(gin_two_down),min_num)],gin_two_down=1),
                                  data.frame(OTU = gin_two_nosig[sample(length(gin_two_nosig),min_num)],gin_two_nosig=1),
                                  data.frame(OTU = gin_two_all[sample(length(gin_two_all),min_num)],gin_two_all=1),
                                  
                                  data.frame(OTU = gin_three_up[sample(length(gin_three_up),min_num)],gin_three_up=1),
                                  data.frame(OTU = gin_three_down[sample(length(gin_three_down),min_num)],gin_three_down=1),
                                  data.frame(OTU = gin_three_nosig[sample(length(gin_three_nosig),min_num)],gin_three_nosig=1),
                                  data.frame(OTU = gin_three_all[sample(length(gin_three_all),min_num)],gin_three_all=1),
                                  
                                  data.frame(OTU = gin_four_up[sample(length(gin_four_up),min_num)],gin_four_up=1),
                                  data.frame(OTU = gin_four_down[sample(length(gin_four_down),min_num)],gin_four_down=1),
                                  data.frame(OTU = gin_four_nosig[sample(length(gin_four_nosig),min_num)],gin_four_nosig=1),
                                  data.frame(OTU = gin_four_all[sample(length(gin_four_all),min_num)],gin_four_all=1),
                                  
                                  data.frame(OTU = qui_two_up[sample(length(qui_two_up),min_num)],qui_two_up=1),
                                  data.frame(OTU = qui_two_down[sample(length(qui_two_down),min_num)],qui_two_down=1),
                                  data.frame(OTU = qui_two_nosig[sample(length(qui_two_nosig),min_num)],qui_two_nosig=1),
                                  data.frame(OTU = qui_two_all[sample(length(qui_two_all),min_num)],qui_two_all=1),
                                  
                                  data.frame(OTU = qui_three_up[sample(length(qui_three_up),min_num)],qui_three_up=1),
                                  data.frame(OTU = qui_three_down[sample(length(qui_three_down),min_num)],qui_three_down=1),
                                  data.frame(OTU = qui_three_nosig[sample(length(qui_three_nosig),min_num)],qui_three_nosig=1),
                                  data.frame(OTU = qui_three_all[sample(length(qui_three_all),min_num)],qui_three_all=1),
                                  
                                  data.frame(OTU = qui_four_up[sample(length(qui_four_up),min_num)],qui_four_up=1),
                                  data.frame(OTU = qui_four_down[sample(length(qui_four_down),min_num)],qui_four_down=1),
                                  data.frame(OTU = qui_four_nosig[sample(length(qui_four_nosig),min_num)],qui_four_nosig=1),
                                  data.frame(OTU = qui_four_all[sample(length(qui_four_all),min_num)],qui_four_all=1),
                                  
                                  data.frame(OTU = not_two_up[sample(length(not_two_up),min_num)],not_two_up=1),
                                  data.frame(OTU = not_two_down[sample(length(not_two_down),min_num)],not_two_down=1),
                                  data.frame(OTU = not_two_nosig[sample(length(not_two_nosig),min_num)],not_two_nosig=1),
                                  data.frame(OTU = not_two_all[sample(length(not_two_all),min_num)],not_two_all=1),
                                  
                                  data.frame(OTU = not_three_up[sample(length(not_three_up),min_num)],not_three_up=1),
                                  data.frame(OTU = not_three_down[sample(length(not_three_down),min_num)],not_three_down=1),
                                  data.frame(OTU = not_three_nosig[sample(length(not_three_nosig),min_num)],not_three_nosig=1),
                                  data.frame(OTU = not_three_all[sample(length(not_three_all),min_num)],not_three_all=1),
                                  
                                  data.frame(OTU = not_four_up[sample(length(not_four_up),min_num)],not_four_up=1),
                                  data.frame(OTU = not_four_down[sample(length(not_four_down),min_num)],not_four_down=1),
                                  data.frame(OTU = not_four_nosig[sample(length(not_four_nosig),min_num)],not_four_nosig=1),
                                  data.frame(OTU = not_four_all[sample(length(not_four_all),min_num)],not_four_all=1)))
    otu[is.na(otu)] <- 0
    rownames(otu) <- otu$OTU
    otu$OTU <- NULL
    phylo_rarefy <- phyloseq(otu_table(otu,taxa_are_rows = T),phy_tree(tree))
    uuf_rarefy <- UniFrac(phylo_rarefy,weighted = F) %>% as.matrix %>% data.frame
    uuf_result[,,i] <- uuf_rarefy %>% as.matrix
    
  }
  rarefy_group <- data.frame(group=c('gin_two_up','gin_two_down','gin_two_nosig','gin_two_all',
                                     'gin_three_up','gin_three_down','gin_three_nosig','gin_three_all',
                                     'gin_four_up','gin_four_down','gin_four_nosig','gin_four_all',
                                     'qui_two_up','qui_two_down','qui_two_nosig','qui_two_all',
                                     'qui_three_up','qui_three_down','qui_three_nosig','qui_three_all',
                                     'qui_four_up','qui_four_down','qui_four_nosig','qui_four_all',
                                     'not_two_up','not_two_down','not_two_nosig','not_two_all',
                                     'not_three_up','not_three_down','not_three_nosig','not_three_all',
                                     'not_four_up','not_four_down','not_four_nosig','not_four_all'))
  rarefy_group <- separate(rarefy_group,col='group',into = c('plant','year','type'),sep='_',remove = F)
  disper <- betadisper(uuf_result[,,1]%>%as.dist,rarefy_group$type)$distances %>% data.frame %>% mutate(.,name=rarefy_group$group,rep=1)
  dis_obj <- betadisper(uuf_result[,,1]%>%as.dist,rarefy_group$type)
  dis_obj_lst <- list()
  dis_obj_lst <- c(dis_obj_lst,list(dis_obj))
  names(disper)[1] <- 'distance'
  for(j in 2:rand){
    disper_test <- betadisper(uuf_result[,,j]%>%as.dist,rarefy_group$type)$distances %>% data.frame %>% mutate(.,name=rarefy_group$group,rep=j)
    dis_obj_test <- betadisper(uuf_result[,,j]%>%as.dist,rarefy_group$type)
    dis_obj_lst <- c(dis_obj_lst,list(dis_obj_test))
    names(disper_test)[1] <- 'distance'
    disper <- rbind(disper,disper_test)
  }
  return(list(uuf_result,disper,dis_obj_lst))
}


my_uuf_resample <- my_rarefy_mean(gin_two_up,gin_two_down,gin_two_nosig,rownames(gin_two),
                                  gin_three_up,gin_three_down,gin_three_nosig,rownames(gin_three),
                                  gin_four_up,gin_four_down,gin_four_nosig,rownames(gin_four),
                                  qui_two_up,qui_two_down,qui_two_nosig,rownames(qui_two),
                                  qui_three_up,qui_three_down,qui_three_nosig,rownames(qui_three),
                                  qui_four_up,qui_four_down,qui_four_nosig,rownames(qui_four),
                                  not_two_up,not_two_down,not_two_nosig,rownames(not_two),
                                  not_three_up,not_three_down,not_three_nosig,rownames(not_three),
                                  not_four_up,not_four_down,not_four_nosig,rownames(not_four),
                                  tree = b_tree,
                                  rand=100)
