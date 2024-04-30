a <- list()
a[[1]] <- c("AATGGTCCTCGATTATTCCCAGGGTGCCGATGAATCTCG")
a[[2]] <- c("AATGGTCCACGATTATTCGCAGGGTGCCGATGAATCTCG")
a[[3]] <- c("AATGGTCCTCGATTATTCCCAGGGTGCCGATGAATCTCG")
a[[4]] <- c("AATGGTCCACGATCATTCCCAGGGTGCAGATGGATCTCG")
a[[5]] <- c("AATGGTCCACGATTATTCGCAGGGTGCCGATGAATCTCG")
a[[6]] <- c("AATGGTCCTCGATTATTCCCAGGGTGCCGATGAATCTCG")
a[[7]] <- c("AATGGTCCTCGATTATTCCCAGGGTGCCGATGAATCTCG")
a[[8]] <- c("AATGGTCCACGATTATTCGCAGGGTGCCGATGAATCTCG")
a[[9]] <- c("AATGGTCCACGATCATTCCCAGGGTGCAGATGGATCTCG")
a[[10]] <- c("AATGGTCCGCGATTATTCTCAGGGTGCGGATGAATCTCG")

a <- lapply(a,function(x){unlist(strsplit(x,split = ""))})
fas_dat <- a
names(fas_dat) <- as.character(c(1:10))


####
###計算Pi
pi_diversity <- c()
step <- if (sliding == 0) window_size else sliding
for (i in seq(1,(2*one_side_length)-window_size+1,by = step)) {
  ##擷取每個區間的長度(1000bp)(temp_list)
  temp_list <- lapply(fas_dat,FUN = sub,pos = i)
  uni <- unique(temp_list)
  ##比對
  freq_list_2 <- list()
  for (p in 1:length(uni)) {
    sub_list <- list()
    freq_list_2[[p]] <- sub_list
  }
  ##將每個個體的名字分配給所屬的不重複序列(freq_list_2) (freq_list為不重複序列本身)
  for (l in 1:10) {
    freq_list <- lapply(uni,FUN = freq,ind = l)
    freq_list_2[[which(freq_list %in% window_size)]] <- append(freq_list_2[[which(freq_list %in% window_size)]],name[l]) 
    
  }
  ##去掉na之後的長度並轉換成比例，用比例當作freq_list_2的名字
  uni_len <- (sapply(freq_list_2,length))/sum(sapply(freq_list_2,length))
  names(freq_list_2) <- uni_len
  
  ##兩兩比較
  pi_value <- list()
  if(length(uni)==1){
    result = 0
  }else{
    for (j in 1:length(uni)) {
      for (k in (j+1):(length(uni))) {
        if(j<length(uni)){
          diff <-which(Map(`!=`,uni[[j]],uni[[k]])==T)
          pi <- length(diff)/length(uni[[j]])
          listname <- paste0(j,",",k)
          pi_value[[listname]] <- pi
          
        }else{
          break
        }
      }
    }
    ###將pi_value list名稱分開，找到各自數字代表的頻率-->freq_value
    list_name <- names(pi_value)
    freq_index <- as.numeric(unlist(sapply(list_name,function(x){strsplit(x,",")})))
    freq_value <- c()
    for (o in 1:length(freq_index)) {
      freq_value[o]<- as.numeric(names(freq_list_2[freq_index[o]]))
    }
    ##將頻率乘上pi值(xixj*pi)
    result <- round(mapply(function(x, y, sublist) x * y * sublist, freq_value[seq(1, length(freq_value), 2)], freq_value[seq(2, length(freq_value), 2)], pi_value),3)
    ###乘上n/n-1
    result <- sum(result)*(length(fas_dat)/(length(fas_dat)-1))
    result <- sum(result)
  }
  pi_diversity <- append(pi_diversity,result)
  pi_value <- lapply(pi_value,function(x){round(pi_value[[1]],2)})