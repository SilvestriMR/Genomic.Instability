MeasureLST <-
function(data, window, ID, workflow=c("SS","MS")){
  if (workflow=="SS"){
    for (x in data[,1]){
      par <- 10000000/window
      ch_name <- x
      idx <- data[data[,1] %in% x,]
      idx1 <- rle(idx[,4])
      idx2 <- grep("-1", idx1$values)
      idx1$lengths[idx2] <- 0 
      c <- as.vector(idx1$lengths)
      c <- append(c, 0)
      c[c<par] <- 0
      c[c>par] <- 1
      count <- as.data.frame(sum(rollsum(c,2) > 1))
      colnames(count) <- x
      assign(x, count)
      
    }

    TabLST <- do.call(cbind, mget(ls(pattern="^chr")))
    TabLST$Score <- rowSums(TabLST)
    rownames(TabLST) <- ID
    TabLST <- TabLST[,c("chr1",
                            "chr2",
                            "chr3",
                            "chr4",
                            "chr5",
                            "chr6",
                            "chr7",
                            "chr8",
                            "chr9",
                            "chr10",
                            "chr11",
                            "chr12",
                            "chr13",
                            "chr14",
                            "chr15",
                            "chr16",
                            "chr17",
                            "chr18",
                            "chr19",
                            "chr20",
                            "chr21",
                            "chr22",
                            "chrX",
                            "Score")]
    TabLST
           
  } else {
      for (y in ID){
        data <- data[[1]]
        wind <- window[1]
        for (x in data[,1]){
          par <- 10000000/wind
          ch_name <- x
          idx <- data[data[,1] %in% x,]
          idx1 <- rle(idx[,4])
          idx2 <- grep("-1", idx1$values)
          idx1$lengths[idx2] <- 0 
          c <- as.vector(idx1$lengths)
          c <- append(c, 0)
          c[c<par] <- 0
          c[c>par] <- 1
          count <- as.data.frame(sum(rollsum(c,2) > 1))
          colnames(count) <- x
          assign(x, count)
          
        }
        
        TabLST <- do.call(cbind, mget(ls(pattern="^chr")))
        TabLST$Score <- rowSums(TabLST)
        rownames(TabLST) <- y
        assign(paste("Bind_",y,sep=""),TabLST)
        data <- data[-1]
        window <- window[-1]
        
      }
    
      Tab_summ <- do.call(rbind, mget(ls(pattern="^Bind_")))
      Tab_summ <- Tab_summ[,c("chr1",
                              "chr2",
                              "chr3",
                              "chr4",
                              "chr5",
                              "chr6",
                              "chr7",
                              "chr8",
                              "chr9",
                              "chr10",
                              "chr11",
                              "chr12",
                              "chr13",
                              "chr14",
                              "chr15",
                              "chr16",
                              "chr17",
                              "chr18",
                              "chr19",
                              "chr20",
                              "chr21",
                              "chr22",
                              "chrX",
                              "Score")]
      Tab_summ
      
    }
  
}
