rm(list=ls(all=TRUE))
draw.letter <- function(letter.vector){
  len <- length(letter.vector)
  for (i in 1:len) {
    cat(ifelse(letter.vector[i] == 1, "x", " ")); # simplify with ifelse() function
    if (i %% 5 == 0) { cat('\n') }
  }
}


mutate<-function(lv, pf){
  pixel<-c(-1,1)
  replace(lv, sample(length(lv),pf, replace= F), sample(pixel,1, replace = F))  
}


signa<-function(x){
  v<-matrix()
  for (i in 1:25){
    if (x[i] > 0){
      v[i] = 1
    }else {v[i] =-1 }
  }
  return(v)
}

D=c(1,1,1,-1,-1,  1,-1,-1,-1,1,  1,-1,-1,-1,1,  1,-1,-1,-1,1,  1,1,1,-1,-1)
J=c(-1,-1,-1,1,1, -1,-1,-1,-1,1,  -1,-1,-1,-1,1,  1,-1,-1,-1,1,  1,1,1,1,1)
C=c(-1,-1,1,1,1,  1,1,-1,-1,-1,  1,1,-1,-1,-1,  1,1,-1,-1,-1, -1,-1,1,1,1 )
M=c(1,-1,-1,-1,1,   1,1,-1,1,1, 1,-1,1,-1,1,  1,-1,-1,-1,1,   1,-1,-1,-1,1)
draw.letter(M)
draw.letter(C)
draw.letter(D)
draw.letter(J)

fv.m=matrix(fv, nrow = 4, byrow = T) 

w=t(fv.m) %*% fv.m
for (i in 1:25) { w[i,i] = 0 }

hopfield<-function(a, mem){
  numm<-length(which(mem != (a)))
  cat("-Mutated letter with",numm, "changed positions-\n")
  draw.letter(a)  
    cat("\nMutations at iterations: ")
  for (r in 1:6){
    a<-w %*% a  
    a<-signa(a)
    a[1]<-1
    ll<-length(which(mem != (a)))
    cat(ll," ")
  }
  cat("\n")
  cat("\n")
  cat("Restored Letter :)\n")
  draw.letter(a)
}
m.M<-mutate(M, 8)
test<-hopfield(m.M,M)



m.J<-mutate(J, 8)
test<-hopfield(m.J,J)

m.C<-mutate(C, 8)
test<-hopfield(m.C,C)

m.D<-mutate(D, 8)
test<-hopfield(m.D,D)




# 
# 
# 
# hopfield<-function(mut, mem){
#   for (r in 1:5){
#     a<-w %*% a  
#     a<-signa(a)
#     ll<-length(which(J != (a)))
#     print(ll)
#   }
#   return(a)
#   #print(signa(out))
#   
# }
# hopfield(d.m,D)
# which(D != d.m)
# 
# activation<-function(x) {
#     av<-vector()
#     for (i in 1:25){
#         if (x[i] != w[i,i]){
#         av[i]<-sum(x[i]* w[i,])
#         } 
#     }
#     return(av)
# }
# 
# activation.5<-function(x){
#   for (r in 1:10){
#   for (i in 1:25){
#       av<-vector()
#       
#       for (j in 1:25){
#       if (x[i] != w[i,j]){
#         av[i]<-sum(x[i]* w[i,])
#       } 
#     }
#   }
#   }
#   return(av)
# }
# 
# signa(activation.5(sadm))
# 
# d.m<-as.matrix(md, byrow= f)
# 
# adm=w %*% d.m
# sadm=signa(adm)
# 
# dim(w)
# 
# which(d.m != sadm)
# draw.letter(sadm)
# draw.letter(d.m)
# M<-mutate(M,4)
# draw.letter(M)
# j.m<-mutate(J,2)
# d.m=mutate(D,8)
# which(D != d.m)
# draw.letter(J)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# which(D != test  )
# igna<-function(x){
#   v<-matrix()
#   for (i in 1:25){
#     if (x[i] > 0){
#       v[i] = 1
#     }else {v[i] =-1 }
#   }
#   return(v)
# }
# 
# 
# 
# hopfield<-function(x){
#   a <- 5
#   for (r in 1:5){
#     #print(ad)
#     a<-2 * a  
#     a<-a+2
#     print(a)
#   }
#   #print(signa(out))
#   
# }
# hopfield(d.m)
# 
# 
# 
#     for (i in 1:25){
#       av<-vector()
#       
#       for (j in 1:25){
#         if (x[j] != w[i,j]){
#           av[i]<-sum(x[j]* w[i,])
#         } 
#       }
#     }
#   }
#   return(av)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# signa.100<-function(x){
#   v<-vector()
#   for (i in 1:100){
#     if (x[i] > 0){
#       v[i] = 1
#     }else {v[i] =-1 }
#   }
#   return(v)
# }
