#############################################################################
#Script to match the taxonomy from the Jetz et al tree to the more updated
#taxonomy from ebird. Code adjusted from online (taxize package) and help from Claude.ai
############################################################################


library(taxize)
library(dplyr)
library(clootl)  # install from: install.packages("clootl", repos = "https://eliotmiller.r-universe.dev")

# ── 1. Load your data ─────────────────────────────────────────────────────────
my_data    <- read.csv("FinalBirds.csv")
jetz_names <- my_data$Species
jetz_names<-gsub("_", " ", jetz_names)

# ── 2. Load the clootl taxonomy and inspect it ────────────────────────────────
# taxonomyGet() is the correct public function — do NOT index clootl_data directly.
# Use the year that matches your target tree (2021–2024 currently supported).
ebird_tax <- taxonomyGet(taxonomy_year = 2024)

# Always inspect before assuming column names
cat("Columns in eBird taxonomy:\n")
print(colnames(ebird_tax))

# ── 3. Resolve your Jetz names against the eBird taxonomy ─────────────────────
# gna_verifier() checks names against global name databases and returns
# the accepted/valid name.
# This makes API calls so may take a while for large datasets.
verified <- gna_verifier(jetz_names)
cat("Columns returned by gna_verifier:\n")
print(colnames(verified))


# The column with the matched accepted name is typically "matched_name"
# but inspect to confirm — it can vary with package version

matched_names <- verified$currentCanonicalSimple

# ── 4. Match resolved names to eBird taxonomy ─────────────────────────────────
# Inspect ebird_tax first to confirm the scientific name column
# It is typically "SCI_NAME" but verify with colnames() above
sci_col <- "SCI_NAME"   # update this after running colnames(ebird_tax)

final_taxonomy_map <- merge(
  x    = data.frame(jetz_name     = verified$submittedName,
                    matched_name  = matched_names),
  y    = ebird_tax,
  by.x = "matched_name",
  by.y = sci_col,
  all.x = TRUE    # keep all your species, even unmatched ones (NAs)
)


# ── 6. Extract the clootl tree for your matched species ───────────────────────
matched_species <- final_taxonomy_map$matched_name[!is.na(final_taxonomy_map$matched_name)]
pruned_tree <- extractTree(
  species       = matched_species,
  label_type    = "scientific",
  taxonomy_year = 2025,
  version = "1.6",
  force = TRUE
)


#Now create database including the species names from the Jetz tree as well as the
#matched names:

names.FinalBirds<-subset(final_taxonomy_map, select=c("matched_name","jetz_name"))
names.FinalBirds[]<-lapply(names.FinalBirds, gsub, pattern=" ", replacement="_")
names.FinalBirds$matched_name[names.FinalBirds$matched_name == ""]<-NA

#Check for duplicated names in matched_name which could result from lumping
anyDuplicated(names.FinalBirds$matched_name[!is.na(names.FinalBirds$matched_name)])

#Extract those duplicates to check them out
dup_vals <- names.FinalBirds$matched_name[
  duplicated(names.FinalBirds$matched_name) | 
    duplicated(names.FinalBirds$matched_name, fromLast = TRUE)
]
dup_vals <- unique(dup_vals[!is.na(dup_vals)])

dup_table_124 <- names.FinalBirds[names.FinalBirds$matched_name %in% dup_vals, ]
dup_table_124 <- dup_table_124[order(dup_table_124$matched_name), ]

dup_table_124

#Drop the duplicated names 
# Identify matched_name values that are duplicated (excluding NA = no match found)
dup_vals <- names.FinalBirds$matched_name[
  duplicated(names.FinalBirds$matched_name) | 
    duplicated(names.FinalBirds$matched_name, fromLast = TRUE)
]
dup_vals <- unique(dup_vals[!is.na(dup_vals)])

# Sanity check: confirm exactly one self-matching row (jetz_name == matched_name) per duplicate group
self_match_counts <- sapply(dup_vals, function(nm) {
  sum(names.FinalBirds$matched_name == nm & names.FinalBirds$jetz_name == nm, na.rm = TRUE)
})
stopifnot(all(self_match_counts == 1))

clean_groups   <- names(self_match_counts[self_match_counts == 1])
problem_groups <- names(self_match_counts[self_match_counts != 1])

rows_to_drop <- (names.FinalBirds$matched_name %in% clean_groups &
                   names.FinalBirds$jetz_name != names.FinalBirds$matched_name) |
  names.FinalBirds$matched_name %in% problem_groups

names.FinalBirds_clean <- names.FinalBirds[!rows_to_drop, ]

anyDuplicated(names.FinalBirds_clean$matched_name[!is.na(names.FinalBirds_clean$matched_name)])  # expect 0

#Final check

nrow(names.FinalBirds) - nrow(names.FinalBirds_clean)  # should equal 14 (dropped synonyms) + however many rows were in that 1 problem group
anyDuplicated(names.FinalBirds_clean$matched_name[!is.na(names.FinalBirds_clean$matched_name)])  # 0


FinalBirds<-merge(my_data, names.FinalBirds_clean, by.x="Species", by.y="jetz_name", all.x=TRUE)
#Should get 0
anyDuplicated(FinalBirds$matched_name[!is.na(FinalBirds$matched_name)])

#Drop rows for which there is no match in matched_name
sum(is.na(FinalBirds$matched_name))
FinalBirds2<-subset(FinalBirds, !is.na(matched_name))

#Select Mass, Lifespan, Migration, Clutch size, Brain size, HWI and Sample size (this is the sample size for Lifespan from AnAge)
FinalBirds2<-subset(FinalBirds, select=c("Species", "matched_name", "Mass", "Migration", "Clutch", "Brain", "HWI", "Lifespan", "Sample.size"))

#Drop NA values in the numeric columns
FinalBirds3<-FinalBirds2[-which(is.na(FinalBirds2[3]) | is.na(FinalBirds2[4]) | is.na(FinalBirds2[5]) | is.na(FinalBirds2[6]) | is.na(FinalBirds2[8])), ]
nrow(FinalBirds3)


#Now log10 transform the numeric variables
#Log10 transform the data
cols<-c("Brain", "Mass", "Lifespan")
FinalBirds3[cols]<-log10(FinalBirds3[cols])
cols2<-c("Migration", "Clutch")
FinalBirds3[cols2]<-log10(FinalBirds3[cols2]+1)

#Transform Sample.size variable to numerical
FinalBirds3$Sample.size<-as.factor(FinalBirds3$Sample.size)
FinalBirds3$Sample.size<-unclass(FinalBirds3$Sample.size)
#Note: given the alphabetical order of the sample size categories (large, medium,
#small, tiny) the numerical order is now 1:large, 2:medium, 3:small, 4:tiny.
names(FinalBirds3)[names(FinalBirds3)=="Sample.size"]<-"Samplesize"

#Load the sample of trees from the two backbones of Jetz et al as well as the 
#sample of trees from Mc Tavish et al. 

JetzHackett<-read.tree("BirdzillaHackett9.tre")
JetzEricson<-read.tree("BirdzillaEricson9.tre")

McTavish<-read.tree("McTavishTreeset.tre")

#Sample 10 trees from the pseudo-posterior distributions of the two Jetz et al backbones
#and from the Mc Tavish tree

set.seed(789)

idx<-sample(length(JetzHackett), 10)  # 10 distinct random positions
JetzHackett10 <- JetzHackett[idx]            

JetzEricson10<-JetzEricson[idx]

idx<-sample(length(McTavish), 10)
McTavish10<-McTavish[idx]

#Compare tip labels of both trees to the two columns of species names

setdiff(FinalBirds3$Species, JetzHackett10[[4]]$tip.label)
#OK
setdiff(FinalBirds3$matched_name, JetzHackett10[[4]]$tip.label)
#15 names differ

setdiff(FinalBirds3$Species, McTavish10[[4]]$tip.label)
#18 names differ

setdiff(FinalBirds3$matched_name, McTavish10[[4]]$tip.label)
#12 names differ, which is strange because these should have been matched above.

#Manual check results:
#Replaced in matched names:
#Ardea ibis by	Bubulcus_ibis
#Astur cooperii by	Accipiter_cooperii
#Astur gentilis	by Accipiter_gentilis
#Botaurus minutus by Ixobrychus_minutus

# Lookup table: old name -> new name
replacements <- c(
  "Ardea_ibis" = "Bubulcus_ibis",
  "Astur_cooperii" = "Accipiter_cooperii",
  "Astur_gentilis" = "Accipiter_gentilis",
  "Botaurus_minutus" = "Ixobrychus_minutus"
)

# Make sure the column is character (a factor would break the assignment)
FinalBirds3$matched_name <- as.character(FinalBirds3$matched_name)

# Replace only the rows that match
hit <- FinalBirds3$matched_name %in% names(replacements)
FinalBirds3$matched_name[hit] <- replacements[FinalBirds3$matched_name[hit]]

#Now need to crop the dataframe to drop the species not present in the McTavish et al tree
setdiff(FinalBirds3$matched_name, McTavish10[[4]]$tip.label)

# Species to drop
drop_species <- setdiff(FinalBirds3$matched_name, McTavish10[[4]]$tip.label)

# Keep only rows whose name is NOT in the list
FinalBirds4 <- FinalBirds3[!(FinalBirds3$matched_name %in% drop_species), ]
#Small sanity check
intersect(FinalBirds4$matched_name, drop_species)  

setdiff(FinalBirds4$matched_name, McTavish10[[4]]$tip.label)

write.csv(FinalBirds4, file="Bird.datMatched8Nov.csv")

#Drop levels for the Samplesize variable (large is 1, medium 2, small 3 and tiny 4)
FinalBirds4$Samplesize <- as.numeric(FinalBirds4$Samplesize)

#Standardize so that the betas of the path coefficients are comparable
FinalBirds4[sapply(FinalBirds4, is.numeric)]<-scale(FinalBirds4[sapply(FinalBirds4, is.numeric)])

write.csv(FinalBirds4, file="Bird.datMatchedScaled.csv")


#Adjust the tip.labels of the Jetz tree to the updated names in matched_names
#Would need first to drop the tips in the tree for species for which there is no match

#Check that trees and dataframe match in the species (matched_name for McTavish tree, and
#Species is for the two Jetz backbones)

tiplabelsTrees<-subset(FinalBirds4, select=c("Species", "matched_name"))

to.drop<-setdiff(JetzHackett10[[9]]$tip.label, tiplabelsTrees$Species)

JetzHackett.pruned<-drop.tip(JetzHackett10, to.drop)
#Check it worked OK
setdiff(JetzHackett.pruned[[9]]$tip.label, tiplabelsTrees$Species)

#Do the same for the Ericson backbone and the Mc Tavish trees

to.drop2<-setdiff(JetzEricson10[[9]]$tip.label, tiplabelsTrees$Species)

JetzEricson.pruned<-drop.tip(JetzEricson10, to.drop2)
#Check it worked OK
setdiff(JetzEricson.pruned[[9]]$tip.label, tiplabelsTrees$Species)

to.drop3<-setdiff(McTavish10[[9]]$tip.label, tiplabelsTrees$matched_name)
McTavish.pruned<-drop.tip(McTavish10, to.drop3)
#Check it worked OK
setdiff(McTavish.pruned[[9]]$tip.label, tiplabelsTrees$matched_name)

#Can now replace the tip.labels in the Jetz backbones with the names of the species
#in the McTavish tree

JetzHacketrenamed<-JetzHackett.pruned
JetzHacketrenamed <- lapply(JetzHacketrenamed, function(tr) {
  tr$tip.label <- tiplabelsTrees$matched_name[match(tr$tip.label, tiplabelsTrees$Species)]
  tr
})
class(JetzHacketrenamed) <- "multiPhylo"

class(JetzHacketrenamed); length(JetzHacketrenamed)          # "multiPhylo", 10
sum(is.na(unlist(lapply(JetzHacketrenamed, `[[`, "tip.label"))))   # 0 = every tip matched
any(sapply(JetzHacketrenamed, function(t) anyDuplicated(t$tip.label) > 0))  # FALSE = no duplicates

setdiff(McTavish.pruned[[9]]$tip.label, JetzHacketrenamed[[9]]$tip.label)

#Repeat for the Ericson backbone

JetzEricsonrenamed<-JetzEricson.pruned
JetzEricsonrenamed<-lapply(JetzEricsonrenamed, function(tr){
  tr$tip.label<-tiplabelsTrees$matched_name[match(tr$tip.label, tiplabelsTrees$Species)]
  tr
})
class(JetzEricsonrenamed)<-"multiPhylo"

setdiff(McTavish.pruned[[9]]$tip.label, JetzEricsonrenamed[[9]]$tip.label)

#Now write the multiphylo objects to a file:

write.tree(McTavish.pruned, file="McTavishmulti.tre")
write.tree(JetzEricsonrenamed, file="JetzEricsonmulti.tre")
write.tree(JetzHacketrenamed, file="JetzHacketmulti.tre")
