### R code from vignette source 'topicmodels.Rnw'

###################################################
### code chunk number 1: topicmodels.Rnw:54-61
###################################################
k <- 30
fold <- 1
options(width = 65, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
library("lattice")
library("topicmodels")
ltheme <- canonical.theme("postscript", FALSE)
lattice.options(default.theme=ltheme)


###################################################
### code chunk number 2: topicmodels.Rnw:510-512
###################################################
cat(paste("R> ", prompt(LDA, filename = NA)$usage[[2]], "\n",
          "R> ", prompt(CTM, filename = NA)$usage[[2]], sep = ""))


###################################################
### code chunk number 3: topicmodels.Rnw:545-552 (eval = FALSE)
###################################################
## control_LDA_VEM <- 
##   list(estimate.alpha = TRUE, alpha = 50/k, estimate.beta = TRUE,
##        verbose = 0, prefix = tempfile(), save = 0, keep = 0,
##        seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
##        var = list(iter.max = 500, tol = 10^-6),
##        em = list(iter.max = 1000, tol = 10^-4),
##        initialize = "random")


###################################################
### code chunk number 4: topicmodels.Rnw:615-621 (eval = FALSE)
###################################################
## control_LDA_Gibbs <- 
##   list(alpha = 50/k, estimate.beta = TRUE,
##        verbose = 0, prefix = tempfile(), save = 0, keep = 0,
##        seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
##        delta = 0.1,
##        iter = 2000, burnin = 0, thin = 2000)


###################################################
### code chunk number 5: topicmodels.Rnw:643-651 (eval = FALSE)
###################################################
## control_CTM_VEM <- 
##   list(estimate.beta = TRUE,
##        verbose = 0, prefix = tempfile(), save = 0, keep = 0,
##        seed = as.integer(Sys.time()), nstart = 1L, best = TRUE,
##        var = list(iter.max = 500, tol = 10^-6),
##        em = list(iter.max = 1000, tol = 10^-4),
##        initialize = "random",
##        cg = list(iter.max = 500, tol = 10^-5))


###################################################
### code chunk number 6: topicmodels.Rnw:720-721
###################################################
data("JSS_papers", package = "topicmodels")


###################################################
### code chunk number 7: topicmodels.Rnw:730-735 (eval = FALSE)
###################################################
## library("OAIHarvester")
## x <- oaih_list_records("https://www.jstatsoft.org/oai", se = "jss:ART")
## x <- x[vapply(x[, "metadata"], length, 1L) > 0L, ]
## JSS_papers  <- oaih_transform(x[, "metadata"])
## JSS_papers <- JSS_papers[order(as.Date(unlist(JSS_papers[, "date"]))), ]


###################################################
### code chunk number 8: topicmodels.Rnw:741-744 (eval = FALSE)
###################################################
## install.packages("corpus.JSS.papers", lib = templib,
##                  repos = "https://datacube.wu.ac.at/", 
##                  type = "source")


###################################################
### code chunk number 9: topicmodels.Rnw:750-753
###################################################
JSS_papers <- JSS_papers[JSS_papers[,"date"] < "2010-08-05",]
JSS_papers <- JSS_papers[sapply(JSS_papers[, "description"], 
                                Encoding) == "unknown",]


###################################################
### code chunk number 10: topicmodels.Rnw:772-774
###################################################
library("tm")
corpus <- Corpus(VectorSource(unlist(JSS_papers[, "description"])))


###################################################
### code chunk number 11: topicmodels.Rnw:783-788
###################################################
Sys.setlocale("LC_COLLATE", "C")
JSS_dtm <- DocumentTermMatrix(corpus, 
   control = list(stemming = TRUE, stopwords = TRUE, minWordLength = 3,
     removeNumbers = TRUE, removePunctuation = TRUE))
dim(JSS_dtm)


###################################################
### code chunk number 12: topicmodels.Rnw:798-807
###################################################
library("slam")
summary(col_sums(JSS_dtm))
term_tfidf <- 
  tapply(JSS_dtm$v/row_sums(JSS_dtm)[JSS_dtm$i], JSS_dtm$j, mean) *
    log2(nDocs(JSS_dtm)/col_sums(JSS_dtm > 0))
summary(term_tfidf)
JSS_dtm <- JSS_dtm[,term_tfidf >= 0.1]
JSS_dtm <- JSS_dtm[row_sums(JSS_dtm) > 0,]
summary(col_sums(JSS_dtm))


###################################################
### code chunk number 13: topicmodels.Rnw:812-813
###################################################
dim(JSS_dtm)


###################################################
### code chunk number 14: topicmodels.Rnw:837-850
###################################################
library("topicmodels")
k <- 30
SEED <- 2010
jss_TM <- 
  list(VEM = LDA(JSS_dtm, k = k, control = list(seed = SEED)),
       VEM_fixed = LDA(JSS_dtm, k = k, 
         control = list(estimate.alpha = FALSE, seed = SEED)),
       Gibbs = LDA(JSS_dtm, k = k, method = "Gibbs",
         control = list(seed = SEED, burnin = 1000, 
           thin = 100, iter = 1000)),
       CTM = CTM(JSS_dtm, k = k, 
         control = list(seed = SEED, 
           var = list(tol = 10^-4), em = list(tol = 10^-3))))


###################################################
### code chunk number 15: topicmodels.Rnw:856-857
###################################################
sapply(jss_TM[1:2], slot, "alpha")


###################################################
### code chunk number 16: topicmodels.Rnw:875-882
###################################################
methods <- c("VEM", "VEM_fixed", "Gibbs", "CTM")
DF <- data.frame(posterior = unlist(lapply(jss_TM, function(x) apply(posterior(x)$topics, 1, max))),
                 method = factor(rep(methods,
                   each = nrow(posterior(jss_TM$VEM)$topics)), methods))
print(histogram(~ posterior | method, data = DF, col = "white", as.table = TRUE,
                xlab = "Probability of assignment to the most likely topic",
                ylab = "Percent of total", layout = c(4, 1)))


###################################################
### code chunk number 17: topicmodels.Rnw:899-902
###################################################
sapply(jss_TM, function(x) 
       mean(apply(posterior(x)$topics, 
                  1, function(z) - sum(z * log(z)))))


###################################################
### code chunk number 18: topicmodels.Rnw:912-913
###################################################
Topic <- topics(jss_TM[["VEM"]], 1)


###################################################
### code chunk number 19: topicmodels.Rnw:917-919
###################################################
Terms <- terms(jss_TM[["VEM"]], 5)
Terms[,1:5]


###################################################
### code chunk number 20: topicmodels.Rnw:928-932
###################################################
(topics_v24 <- 
 topics(jss_TM[["VEM"]])[grep("v024", vapply(JSS_papers[, "identifier"], 
                                             "[", 2, FUN.VALUE = ""))])
most_frequent_v24 <- which.max(tabulate(topics_v24))


###################################################
### code chunk number 21: topicmodels.Rnw:939-940
###################################################
terms(jss_TM[["VEM"]], 10)[, most_frequent_v24]


