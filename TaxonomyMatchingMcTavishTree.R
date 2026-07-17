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
dup_vals <- names.FinalBirds_intree$matched_name[
  duplicated(names.FinalBirds_intree$matched_name) | 
    duplicated(names.FinalBirds_intree$matched_name, fromLast = TRUE)
]
dup_vals <- unique(dup_vals[!is.na(dup_vals)])

dup_table_124 <- names.FinalBirds_intree[names.FinalBirds_intree$matched_name %in% dup_vals, ]
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

#Select only Mass, Lifespan, Migration, Clutch size and Brain size
FinalBirds2<-subset(FinalBirds, select=c("Species", "matched_name", "Mass", "Lifespan", "Migration", "Clutch", "Brain"))

#Add tarsus length and hand wing index to the database (23 June 2026)

avonet<-read.csv("AVONETBirdtree.csv")
avonet<-subset(avonet, select=c("Species3", "Hand.Wing.Index"))
FinalBirds3<-merge(FinalBirds2, avonet, by.x="Species", by.y="Species3", all.x = TRUE)


#Drop NA values in the numeric columns
FinalBirds3<-FinalBirds2[-which(is.na(FinalBirds2[3]) | is.na(FinalBirds2[4]) | is.na(FinalBirds2[5]) | is.na(FinalBirds2[6]) | is.na(FinalBirds2[7])), ]
nrow(FinalBirds3)


#Now log10 transform the numeric variables
#Log10 transform the data
cols<-c("Brain", "Mass", "Lifespan")
FinalBirds3[cols]<-log10(FinalBirds3[cols])
cols2<-c("Migration", "Clutch")
FinalBirds3[cols2]<-log10(FinalBirds3[cols2]+1)

#Standardize so that the betas of the path coefficients are comparable
FinalBirds3[sapply(FinalBirds3, is.numeric)]<-scale(FinalBirds3[sapply(FinalBirds3, is.numeric)])

write.csv(FinalBirds3, file="Bird.datMatchedScaled.csv")

#Adjust the tip.labels of the Jetz tree to the updated names in matched_names
#Would need first to drop the tips in the tree for species fro which there is no match

tiplabelsTrees<-subset(FinalBirds2, select=c("Species", "matched_name"))
tiplabelsTrees<-na.omit(tiplabelsTrees)

#Prune the consensus Jetz tree to include only the species in this pruned vector
#can use this to replace the species names in the tree. Then I can prune the new
#tree to match the number of species in the Mc Tavish tree I downloaded
library(phangorn)

Jetztrees<-read.tree("BirdzillaHackett9.tre")
Jetzmaxcred<-maxCladeCred(Jetztrees)

setdiff(Jetzmaxcred$tip.label, tiplabelsTrees$Species)
#There are 123 species in the tree without matching species in the dataframe

to.drop<-setdiff(Jetzmaxcred$tip.label, tiplabelsTrees$Species)

Jetzmc.pruned<-drop.tip(Jetzmaxcred, to.drop)

#Can now replace the tip.labels in the Jetzmaxcred with the names of the species
#in the McTavish tree

Jetzrenamed<-Jetzmc.pruned
Jetzrenamed$tip.label <- tiplabelsTrees$matched_name[match(Jetzmc.pruned$tip.label,tiplabelsTrees$Species)]

#Check whether both trees match
mctav_phy<-read.tree("McTavishtree.tre")

setdiff(mctav_phy$tip.label, Jetzrenamed$tip.label)
#need to drop "Setophaga_pinus" from McTavish

mctav_phy<-drop.tip(mctav_phy, "Setophaga_pinus")

to.drop<-setdiff(Jetzrenamed$tip.label, mctav_phy$tip.label)
#The 623 species that did not match between both trees from the start

Jetzrenamed<-drop.tip(Jetzrenamed, to.drop)

setdiff(Jetzrenamed$tip.label, mctav_phy$tip.label)

write.tree(Jetzrenamed, file="JetzMCC_matched.tre")

write.tree(mctav_phy, file="McTavishtree.tre")

#Now make a multiphylo object with both trees
multi.phyNew<-c(Jetzrenamed, mctav_phy)

write.nexus(multi.phyNew, file="MultiPhyJetzMcTav.nex")

#Need to make sure the names in the datafile with data for all species match 
#the names in the new multiphylo tree

bird.dat<-read.csv("Bird.datMatchedScaled.csv")
bird.phy.cut<-read.nexus("MultiPhyJetzMcTav.nex")

setdiff(bird.dat$matched_name, bird.phy.cut[[1]]$tip.label)
#Need to drop these 15 species from the dataframe

cut.data<-setdiff(bird.dat$matched_name, bird.phy.cut[[2]]$tip.label)
bird.dat2<-subset(bird.dat, !(matched_name %in% cut.data))

setdiff(bird.dat2$matched_name, bird.phy.cut[[1]]$tip.label)

write.csv(bird.dat2, file="Bird.datMatchedScaled.csv")
