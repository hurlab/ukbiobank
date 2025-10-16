# Filter on diabetes


df$diabetes=0
df$diabetes[df$f.120007.0.0==1]=1
df_diabetes=df[df$diabetes==1,]

# Drop questionnaire rows lacking an ID
if ("ID" %in% names(df_diabetes)) {
  df_diabetes=df_diabetes[!is.na(df_diabetes$ID), ]
}

# Clean MNSI
df_diabetes$f.120071.0.0[df_diabetes$f.120071.0.0<0]=NA
df_diabetes$f.120072.0.0[df_diabetes$f.120072.0.0<0]=NA
df_diabetes$f.120073.0.0[df_diabetes$f.120073.0.0<0]=NA
df_diabetes$f.120074.0.0[df_diabetes$f.120074.0.0<0]=NA
df_diabetes$f.120075.0.0[df_diabetes$f.120075.0.0<0]=NA
df_diabetes$f.120076.0.0[df_diabetes$f.120076.0.0<0]=NA
df_diabetes$f.120077.0.0[df_diabetes$f.120077.0.0<0]=NA
df_diabetes$f.120078.0.0[df_diabetes$f.120078.0.0<0]=NA
df_diabetes$f.120079.0.0[df_diabetes$f.120079.0.0<0]=NA
df_diabetes$f.120080.0.0[df_diabetes$f.120080.0.0<0]=NA
df_diabetes$f.120081.0.0[df_diabetes$f.120081.0.0<0]=NA
df_diabetes$f.120082.0.0[df_diabetes$f.120082.0.0<0]=NA
df_diabetes$f.120083.0.0[df_diabetes$f.120083.0.0<0]=NA
df_diabetes$f.120084.0.0[df_diabetes$f.120084.0.0<0]=NA
df_diabetes$f.120085.0.0[df_diabetes$f.120085.0.0<0]=NA

# Create 2 questions that are "reversed" in the MNSI 
df_diabetes$f.120083.0.0_rev=NA
df_diabetes$f.120083.0.0_rev[df_diabetes$f.120083.0.0==1]=0
df_diabetes$f.120083.0.0_rev[df_diabetes$f.120083.0.0==0]=1


df_diabetes$f.120077.0.0_rev=NA
df_diabetes$f.120077.0.0_rev[df_diabetes$f.120077.0.0==1]=0
df_diabetes$f.120077.0.0_rev[df_diabetes$f.120077.0.0==0]=1

# Calculate the MNSI Score
df_diabetes$MNSI=df_diabetes$f.120071.0.0+
  df_diabetes$f.120072.0.0+
  df_diabetes$f.120073.0.0+
  # Q4 Not scored
  df_diabetes$f.120075.0.0+
  df_diabetes$f.120076.0.0+
  # Q7 Needs to be reversed
  df_diabetes$f.120077.0.0_rev+
  df_diabetes$f.120078.0.0+
  df_diabetes$f.120079.0.0+
  # Q10 Not scored
  df_diabetes$f.120081.0.0+
  df_diabetes$f.120082.0.0+
  # Q13 Needs to be reversed
  df_diabetes$f.120083.0.0_rev+
  df_diabetes$f.120084.0.0+
  df_diabetes$f.120085.0.0

# Allow partial sums to stand when they already exceed the PN threshold
mnsi_items=c(
  "f.120071.0.0",
  "f.120072.0.0",
  "f.120073.0.0",
  "f.120075.0.0",
  "f.120076.0.0",
  "f.120077.0.0_rev",
  "f.120078.0.0",
  "f.120079.0.0",
  "f.120081.0.0",
  "f.120082.0.0",
  "f.120083.0.0_rev",
  "f.120084.0.0",
  "f.120085.0.0"
)
mnsi_partial=rowSums(df_diabetes[, mnsi_items, drop=FALSE], na.rm=TRUE)
mnsi_non_missing=rowSums(!is.na(df_diabetes[, mnsi_items, drop=FALSE]))
mnsi_can_impute=is.na(df_diabetes$MNSI) & mnsi_partial>=3 & mnsi_non_missing>0
df_diabetes$MNSI[mnsi_can_impute]=mnsi_partial[mnsi_can_impute]

# Create PN Outcome
df_diabetes$pn=NA
df_diabetes$pn[df_diabetes$MNSI>3]=1
df_diabetes$pn[df_diabetes$MNSI<=3]=0

