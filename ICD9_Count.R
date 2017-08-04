#Computes the average ICD9s per HL / average ICD9s per NH

HL_Only <- ICD9.df[ICD9.df$PID %in% cases$PID,]
length(unique(HL_Only$PID))
NH_Only <- ICD9.df[ICD9.df$PID %in% controls$PID,]
length(unique(NH_Only$PID))

HL_ICD9_Count <- sum(HL_Only$Code_Count)
NH_ICD9_Count <- sum(NH_Only$Code_Count)

(HL_ICD9_Count/length(unique(HL_Only$PID)))/(NH_ICD9_Count/length(unique(NH_Only$PID)))

mean(aggregate(Code_Count ~ PID, data = HL_Only,sum)$Code_Count)/mean(aggregate(Code_Count ~ PID, data = NH_Only,sum)$Code_Count)

HL_CountPerPID <- aggregate(Code_Count ~ PID, data = HL_Only,sum)$Code_Count
NH_CountPerPID <- aggregate(Code_Count ~ PID, data = NH_Only,sum)$Code_Count

t.test(HL_CountPerPID,NH_CountPerPID)$p.value
