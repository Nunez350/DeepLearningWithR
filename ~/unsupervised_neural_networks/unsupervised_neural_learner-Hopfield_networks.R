rm(list=ls(all=TRUE))
#draws pixalated letter
draw.letter <- function(letter.vector){
  len <- length(letter.vector)
  for (i in 1:len) {
    cat(ifelse(letter.vector[i] == 1, "x", " ")); # simplify with ifelse() function
    if (i %% 5 == 0) { cat('\n') }
  }
}

#mutates desired number of letter pixels
mutate<-function(lv, pf){
  pixel<-c(-1,1)
  replace(lv, sample(length(lv),pf, replace= F), sample(pixel,1, replace = F))  
}

#activation dynamics
#returns signa value from activity values
signa<-function(x){
  v<-matrix()
  for (i in 1:25){
    if (x[i] > 0){
      v[i] = 1
    }else {v[i] =-1 }
  }
  return(v)
}
#applies hopfield network to compare synaptic connections and correct mutated pixels
hopfield<-function(a, mem){
  numm<-length(which(mem != (a)))
  cat("-Mutated letter with",numm, "changed positions-\n")
  draw.letter(a)  
  cat("\nMutations at iterations: ")
  for (r in 1:15){
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

D=c(1,1,1,1,-1,  -1,1,-1,-1,1,  -1,1,-1,-1,1,  -1,1,-1,-1,1,  -1,1,1,1,-1)
J=c(1,1,1,1,1, -1,-1,-1,1,-1,  -1,-1,-1,1,-1,  1,-1,-1,1,-1,  1,1,1,-1,-1)
C=c(-1,1,1,1,1,  1,-1,-1,-1,-1,  1,-1,-1,-1,-1,  1,-1,-1,-1,-1, -1,1,1,1,1 )
M=c(1,-1,-1,-1,1,   1,1,-1,1,1, 1,-1,1,-1,1,  1,-1,-1,-1,1,   1,-1,-1,-1,1)

draw.letter(M)
draw.letter(C)
draw.letter(D)
draw.letter(J)
#make I and t


#flattens and combines all letter vectors
fv<-c(D,J,C,M)
fv.m=matrix(fv, nrow = 4, byrow = T) 
#Hebbian weights change by inner product of activation states
#inner product/dot product function to calculate weight matrix
w=t(fv.m) %*% fv.m
#eliminates activation states for neurons that are connected to theselves
for (i in 1:25) { w[i,i] = 0 }


m.M<-mutate(M, 5)
test<-hopfield(m.M,M)


m.J<-mutate(J, 5)
test<-hopfield(m.J,J)

m.C<-mutate(C, 5)
test<-hopfield(m.C,C)

m.D<-mutate(D, 10)
test<-hopfield(m.D,D)








