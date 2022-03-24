# This script is to subset the full EntoGEM library by taxonomic tags to only
# the odonata articles, and to then fit topic models based on the results of 
# Phase 1 screening to predict probability of inclusion for each study.

# Note that the datasets needed to run this script are NOT included. The underlying 
# data are copyrighted abstracts which we cannot archive. Please contact the 
# lead and/or corresponding authors (Eliza Grames and Chris Elphick) for data.

# Functions and libraries ------------------------------------------------------
library(dplyr) # need to load for pipes

# This function pulls out just good keywords from the abstracts
remove_junk <- function(x){
  tmp <- litsearchr::extract_terms(x, method = "fakerake", min_freq = 1, 
                                   min_n = 1, ngrams = FALSE)
  tmp <- gsub(" ", "_", tmp)
  return(paste(tmp, collapse = " "))
}

# Read in the results of Phase 1 screening -------------------------------------

articles <- read.csv("./Articles_P16612_0323_A7091.csv")
answers <- read.csv("./Answers_P16612_0323_A7091.csv")
users <- read.csv("./UserAnswers_P16612_0323_A7091.csv")

# There are 15 articles which are still flagged as a conflict
# These were third-party screened by GAM, but for some reason were not
# coded as resolutions, so we need to resolve those manually
conflicts <- answers$Article.ID[answers$Status=="conflict"]

for(i in 1:length(conflicts)){
  if(i %in% c(1:8,10,11,13:15)){
    resolutions <- users$Include[(users$Article.ID==conflicts[i] & users$User.Name=="montgomery.graham")]
    answers$Include[answers$Article.ID==conflicts[i]] <- resolutions
    rm(resolutions)
  }    
}

# There are still two with missing resolutions; default to inclusion
answers$Include[answers$Include=="false;;;true"] <- "true"

# Now we can merge the abstracts to the user decisions
# We need the terms in the abstracts for the topic models
answers <- dplyr::left_join(answers, articles[,c('Article.ID', 'Abstract')], 
                            by="Article.ID")

# Read in EntoGEM library -----------------------------------------------------

# Note that this version is actually only partially deduplicated.
# This is because EMG did the original taxonomic tagging before running
# a more sophisticated deduping algorithm. It takes quite a while to run, so
# for consistency, we now take taxonomic subsets from the partially deduped 
# EntoGEM library, and then deduplicate those.
load("./partial_dedupe.rda")
doi_dedupe$title <- gsub(" \\[electronic resource\\]", "", doi_dedupe$title)

# Read in the taxa tags --------------------------------------------------------

load("./EntoGEM_insect_dfm.rda")
presence_only <- output # we only care if an order is present or absent
presence_only[presence_only>0] <- 1
orders <- colSums(presence_only)
barplot(sort(orders), horiz=T, las=1) # lots of leps

# For the taxa-specific projects, we only include studies that match *only*
# to that taxa; others are shuffled into community projects
multiples <- output[rowSums(presence_only)>1,]
single_taxa <- output[rowSums(presence_only)==1,]
single_taxa <- data.frame(single_taxa)

odonates <- gsub("text", "", rownames(single_taxa[single_taxa$Odonata>0,]))

dragons <- doi_dedupe[as.numeric(odonates),]

# Remove articles with exact title and abstract matches when NA is not considered
exact_title <-
  which(duplicated(tolower(trimws(
    tm::removePunctuation(dragons$title)
  ))) & !is.na(dragons$title))
exact_abstract <-
  which(duplicated(tolower(trimws(
    tm::removePunctuation(dragons$abstract)
  ))) & !is.na(dragons$abstract))

dragons <- dragons[-unique(append(exact_title, exact_abstract)), ]

# Prep text for models ---------------------------------------------------------

new_data <- data.frame(Article.ID=dragons$id, Status=rep(NA, nrow(dragons)), 
                       Include=rep(NA, nrow(dragons)),
                       Title=dragons$title, Abstract=dragons$abstract)

answers <- answers[,colnames(answers)%in%c(colnames(new_data))]

modeldat <- rbind(answers, new_data)
modeldat$text <- tm::removeNumbers(tolower(paste(modeldat$Title, modeldat$Abstract)))
modeldat$text <- gsub("species", "", modeldat$text) # swamps models

# Pull text fields from bib data
tiabs <- tolower(modeldat$text)

# Remove irrelevant terms, stopwords, and publishing terms
tiabs <- unlist(lapply(tiabs, remove_junk))
# save(tiabs, file="./cleaned_tiabs.rda") # save this because it takes a while

modeldat$text <- tiabs

modeldat$text <- gsub("abstract_authors", "", modeldat$text) # not really a term

# Assign article classifications -----------------------------------------------

# Sort screening decisions into yes, no, and maybe
modeldat$classification <- rep("maybe", nrow(modeldat)) # default to maybe

# We only want to train on the articles that definitely meet criteria
modeldat$classification[modeldat$Include=="true"] <- "yes"
modeldat$classification[grep("unclear", modeldat$Multi.year.study.)] <- "maybe"
modeldat$classification[modeldat$Include=="false"] <- "no" # Must run last
modeldat$classification[(nrow(answers)+1):nrow(modeldat)] <- NA

modeldat$docID <- as.numeric(rownames(modeldat))

doctext <- dplyr::tibble(modeldat[,c('text','classification','docID')])

tidy_doctext <- doctext %>%
  tidytext::unnest_tokens(word, text) %>%
  group_by(word) %>%
  filter(n() > 10) %>%
  ungroup()

doctext_split <- tidy_doctext %>%
  select(docID) %>%
  rsample::initial_split()

# There are no classifications for the odonate articles
train_data <- tibble(data.frame(docID=modeldat$docID[!is.na(modeldat$classification)]))
new_data <- tibble(data.frame(docID=modeldat$docID[is.na(modeldat$classification)]))

sparse_words <- tidy_doctext %>%
  count(docID, word) %>%
  inner_join(train_data) %>%
  tidytext::cast_sparse(docID, word, n)

word_rownames <- as.integer(rownames(sparse_words))

doctext_joined <- data_frame(docID = word_rownames) %>%
  left_join(doctext %>%
              select(docID, classification))

# Run topic models -------------------------------------------------------------
doMC::registerDoMC(cores = 8)

is_included <- doctext_joined$classification == "yes"

model <- glmnet::cv.glmnet(sparse_words, is_included,
                   family = "binomial",
                   parallel = TRUE, keep = TRUE)

# Make predictions for the odonata articles ------------------------------------

coefs <- model$glmnet.fit %>%
  broom::tidy() %>%
  filter(lambda == model$lambda.1se)

intercept <- coefs %>%
  filter(term == "(Intercept)") %>%
  pull(estimate)

classifications <- tidy_doctext %>%
  inner_join(new_data) %>%
  inner_join(coefs, by = c("word" = "term")) %>%
  group_by(docID) %>%
  summarize(score = sum(estimate)) %>%
  mutate(probability = plogis(intercept + score))

hist(classifications$probability, 100)
sum(classifications$probability>.2)

# Note that these predictions are different from the original model run for
# the odonates, which we ran prior to completing Phase 1 so there was a 
# different training set. The models and predictions evolve as the project
# progresses. 
priorities <- classifications$docID[classifications$probability>0.2]
reserve <- classifications$document[classifications$probability<=0.2]

# Write results as bibliographic data ------------------------------------------
sendtoris <- dragons[which((dragons$title %in% modeldat$Title[priorities])),]

colnames(sendtoris)[8] <- "author"
sendtoris$type <- rep("JOUR", nrow(sendtoris))

bibdat <- synthesisr::as.bibliography(sendtoris)
bibdat <- synthesisr::write_ris(bibdat)

