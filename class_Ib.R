getwd()
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")


# Load the dataset from the raw_data folder
patient_info <- read.csv("raw_data/patient_info.csv")

# Inspect the structure of the dataset
str(patient_info)
summary(patient_info)
head(patient_info)
tail(patient_info)

# Convert variables to appropriate data types
# 'gender', 'diagnosis', and 'smoker' are categorical and should be factors
patient_info$gender <- as.factor(patient_info$gender)
patient_info$diagnosis <- as.factor(patient_info$diagnosis)
patient_info$smoker <- as.factor(patient_info$smoker)

# Create a new binary factor for smoking status (1 for "Yes", 0 for "No")
# We will create a new column 'smoker_binary' where 1 = "Yes" and 0 = "No".
patient_info$smoker_binary <- ifelse(patient_info$smoker == "Yes", 1, 0)
patient_info$smoker_binary <- as.factor(patient_info$smoker_binary)

# Inspect the cleaned data structure to verify changes
cat("\n--- Cleaned Data Structure ---\n")
str(patient_info)

# Save the cleaned dataset
# To be saved in the "clean_data" folder.
write.csv(patient_info, file = "clean_data/patient_info_clean.csv", row.names = FALSE)
