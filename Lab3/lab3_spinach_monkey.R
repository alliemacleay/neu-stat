
library(expm)

print("Finding MLE.Score and Distance")

f <- read.table("lab3_f.txt")
pam <- read.table("macleay_lab3_pam.txt")

# just to initialize d as an object
d <- read.table("macleay_lab3_pam.txt")

# divide pam matrix by 10, 000
pam1 <- pam / 10000
# Spinach sequence
Seq1 <- 'NGTKESITKLVSDLNSATLEADVDVVVAPPFVYIDQVKSSLTGRVEISAQNCW'
# Monkey sequence
Seq2<-'NGRKQNLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCY'
# 20 Amino Acids
AA <- 'ARNDCQEGHILKMFPSTWYV'

# initilize S1 and S2
S1 <- c() # Made a vector so indexing below is correct. For a matrix, you
S2 <- c() # always need a comma in the brackets, or you can get unexpected
          # results.
          # A matrix index : M[i, ]
          # A vector index : M[i]

# Translate each sequence into an array of number values
for(i in 1:nchar(Seq1)) {
  ch <- substring(Seq1, i, i)
  S1[i] <- regexpr(ch, AA)[1]
}

for(i in 1:nchar(Seq2)) {
  ch <- substring(Seq2, i, i)
  S2[i] <- regexpr(ch, AA)[1]
}

# initialize ProbSum
ProbSum <- c() # Made a vector

for(k in 1:150) {
  # pam matrix to the power of k
  pamk <- as.matrix(pam1) %^% k # Need to do matrix mult here

  for(i in 1:20) {
    for(j in 1:20) {
      # populate d with the pam value divided by
      # the frequency
      d[i, j] <- pamk[i, j] / f[i, 1]
      # Is this a mistake in the lab equation?
      # below is the equation I have from my notes
      # last fall.
      # log( P[i, j] / (Qi * Qj) )
      # - in R -
      # d[i, j] <- pamk[i, j]/(f[i, 1]*f[j, 1])
    }
  }

  # * Normalize the values *
  # The following line is where I have a problem
  # log10 is returning negative values
  D <- 10*log10((d+t(d))/2)
  # matrix + the transform of the matrix divided by 2
  # normalizes the matrix to ensure it is symmetrical
  m <- as.matrix(D)

  # * Remove -Inf values returned for log10 of 0
  m[!is.finite(m)] <- 0
  D <- as.data.frame(m)
  ProbSum[k] <- 0
  for(i in 1:length(S1)) {
    # * for each amino acid return the normalized score
    # of transitioning from one to the next from the matrix
    # and sum all values
    ProbSum[k] <- D[S1[i], S2[i]] + ProbSum[k]
  }
}

# * find the max likelihood estimation
MLE.score <- max(ProbSum)

# * find the PAM distance from the max score
MLE.PAM.distance <- which(ProbSum==MLE.score)

# * set up x and y to graph values
x <- 15 * (1:10) # I just vectorized this, more efficient than for-loop
y <- ProbSum[x]

print("MLE.score and MLE.PAM.distance are populated.  Run plot(x, y).")

plot(x, y)
