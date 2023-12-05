
setwd("C:/Users/HP OMEN/Desktop/Drmac")
dico <- unname(as.matrix(read.table("MYdico.txt"))) #matrica bridova
index <- read.table("MYindex.txt")$V1 #rije?i poredane po indeksu
uniqueVertices <- sort(union(dico[,1],dico[,2])) #skup indeksa 
nverts <- length(uniqueVertices) #ukupan broj vrhova
nedges <- nrow(dico) #ukupan broj bridova

#frekvencije pojavljivanja svake rije?i u grafu
uniqueVertices <- cbind(uniqueVertices, matrix(0, nverts,3))
dim(uniqueVertices)

for(i in 1:nedges){ #petlji treba 20ak sekundi na mom laptopu
  i1 <- dico[i,1]
  i2 <- dico[i,2]
  uniqueVertices[i1,2] <- uniqueVertices[i1,2]+1 #drugi stupac - koliko je rije? puta pokazivala na drugu
  uniqueVertices[i2,3] <- uniqueVertices[i2,3]+1 #tre?i stupac - koliko je rije? puta bila pokazana od druge
  uniqueVertices[c(i1,i2),4] <- uniqueVertices[c(i1,i2),4]+1 #?etvrti stupac - njihova suma
}

cbind(index[order(uniqueVertices[,4],decreasing = T)[1:185]]) #lista 50 "najpopularnijih" rije?i - ve?ina veznici, ?estice,...
rownames(uniqueVertices) <- index #imenuj vrhove
#neki vrhovi se pojavljuju pre ?esto - veznici, ?estice. Stavit ?emo da takve rije?i nisu nikome sinonimi
Banned <- uniqueVertices[(uniqueVertices[,4]) > 784,1] 
length(Banned) #ima ih 183

# R nema implementiran binary search - ubrzava algoritam dosta
binary_search <- function(val, vertices){
  r <- length(vertices)
  l <- 1
  while(r >= l){
    mid = (r+l) %/% 2
    if(vertices[mid] == val){return(mid)}
    else if(val > vertices[mid]){l=mid+1}
    else{r = mid-1}
  }
  return(NA)
}

#funkcija koja ra?una B%*%Z%*%t(A)+t(B)%*%Z%*%A
#ideja je da ne moramo ra?unati matricu susjedstva B ako je pograf velik (samo podgraf) koji je dan preko bridova
#stavio sam da ne ra?unamo B osim ako je podgraf ima vi?e od 5000 vrhova (matrica sa 25e6 elemenata)
Do_iteration <- function(inds1,inds2,Z){ #prima indekse podgrafa i Z
  n <- nrow(Z)
  Z1 <- cbind(Z[,2], Z[,3],rep(0,n)) #matrica Z%*%t(A)
  Z2 <- cbind(rep(0,n),Z[,1],Z[,2]) #matrica Z%*%A
  nedges_small <- length(inds1) #broj bridova podgrafa
  Z10 <- matrix(0,n,3)
  Z20 <- matrix(0,n,3)
  for(i in 1:nedges_small){
    Z10[inds1[i],] <- Z10[inds1[i],] + Z1[inds2[i],]
    Z20[inds2[i],] <- Z20[inds2[i],] + Z2[inds1[i],]
  }
  Z <- Z10 + Z20
  Z1 <- cbind(Z[,2], Z[,3],rep(0,n)) #matrica Z%*%t(A)
  Z2 <- cbind(rep(0,n),Z[,1],Z[,2]) #matrica Z%*%A
  Z10 <- matrix(0,n,3)
  Z20 <- matrix(0,n,3)
  for(i in 1:nedges_small){
    Z10[inds1[i],] <- Z10[inds1[i],] + Z1[inds2[i],]
    Z20[inds2[i],] <- Z20[inds2[i],] + Z2[inds1[i],]
  }
  return(Z10 + Z20)
}


#Sama funkcija
Synonym_reader <- function(word, n_list = 10, tol = 5*.Machine$double.eps){ #izbacuje rang listu duljine n_list
  t1 <- Sys.time()#po?ni mjeriti vrijeme
  w_ind <- which(index == word) #indeks rije?i
  ins <- dico[dico[,1] == w_ind, 2] #indeksi rije?i na koje word pokazuje
  outs <- dico[dico[,2] == w_ind, 1] #indeksi rije?i koji pokazuju na word
  vertices <- sort(setdiff(union(outs,ins),Banned))#indeksi vrhova podgrafa - izbacuje duplikate
  n <- length(vertices) #broj vrhova u podgrafu
  print(paste("Povezan sa",n,"vrhova"))
  if(n==0){#samo ako rije? nije u dico - vra?a NULL
    return(NULL)  
  }
  wor_vertices <- index[vertices] #imena vrhova grafa
  if(length(vertices) == 1){ #ako imamo samo jedan vrh koji nije tra?ena rije? - npr za "phone"
    return (unname(cbind(wor_vertices[1])))
  }
  dico_small <- dico[dico[,1]%in%vertices & dico[,2]%in%vertices,] #bridovi u podgrafu
  nedges_small <- nrow(dico_small) #broj bridova u podgrafu
  print(paste("Broj bridova podgrafa:",nedges_small))
  if(nedges_small == 0){#ako je podgraf prazan, vra?amo null - rije?i s kojima je povezan, me?usobno nisu povezane
    return(NULL)  
  }
  
  #mmm <- matrices(dico_small)
  
  inds1 <- sapply(dico_small[,1],function(x) binary_search(x,vertices))
  inds2 <- sapply(dico_small[,2],function(x) binary_search(x,vertices))
  Z <- matrix(1,n,3)
  #Z <- cbind(rep(1, n) ,rep(1, n) ,rep(1, n) ) #matrica jedinica
  Ztemp <- Z #pomo?na matrica - za kriterij zaustavljanja
  psi <- Inf #za kriterij zaustavljanja
  print("Zapo?ni iteracije")
  j = 0 #broja? iteracija
  while(T){ #algoritam - pametan dio koda
    j <- j+1
    Z <- Do_iteration(inds1,inds2,Z)
    Z <- Z / norm(Z, type = "F")
    psi <- sum((Z-Ztemp)^2)
    if(psi < tol){
      break  
    }
    Ztemp <- Z
  }
  print(paste("Ukupno",j,"iteracija"))
  row.names(Z) <- wor_vertices #nazovi retke po rije?ima
  Z <- cbind(Z[order(Z[,2], decreasing = T),]) #sortiraj po drugom stupcu
  rijeci <- row.names(Z)
  t2 <- Sys.time() #zavr?i mjerenje
  print(t2-t1)
  return(cbind(rijeci[1:min(n,n_list)]))
}



#upi?i rije?
Synonym_reader("disappear")
Synonym_reader("bird")
Synonym_reader("sleep")
Synonym_reader("meal")
Synonym_reader("rectangle")



