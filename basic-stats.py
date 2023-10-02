import pandas as pd
import sys
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from scipy import stats
from statistics import mean 
# Make sure you have the above libraries 

# RUN SCRIPT IN TERMINAL BY ENTERING "python basic-stats.py data.csv":
data = pd.read_csv(sys.argv[1], header=None, names=['Bind', 'Peptide'])

# Display the first few rows of the data
# print(data.head())

# Check the distribution of positive and negative examples
label_counts = data['Bind'].value_counts()
print("Positive examples:", label_counts[1])
print("Negative examples:", label_counts[0])
print("Peptide length statistics:")
data['Peptide Length'] = data['Peptide'].apply(len)

# Add peptide lengths column
print(data['Peptide Length'].describe())

amino_acid_list = 'ACDEFGHIKLMNPQRSTVWY'

# === ONE-HOT ENCODING ANALYSIS ===
# Define a function to one-hot encode amino acids
def one_hot_encoding(sequence):
    encoding = [1 if aa in sequence else 0 for aa in amino_acid_list]
    return encoding

# Add one-hot encoded column table for amino acid composition
amino_acid_features = data['Peptide'].apply(one_hot_encoding).tolist()
amino_acid_df = pd.DataFrame(amino_acid_features, columns=list(amino_acid_list))
data = pd.concat([data, amino_acid_df], axis=1)

# Display the modified dataset with amino acid composition features
# print(data.head())

# Perform a logistic regression analysis
# Split the data into training and testing sets
X = data.drop(['Bind', 'Peptide', 'Peptide Length'], axis=1)  # Use one-hot encoded table as features
y = data['Bind']  # Target variable

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a logistic regression classifier
clf = LogisticRegression(random_state=42)
clf.fit(X_train, y_train)

# Make predictions on the test set
y_pred = clf.predict(X_test)

# Evaluate the classifier
accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)
print("Classification Report:\n", classification_report(y_test, y_pred))
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

# === HYDROPHOBICITY ANALYSIS ===
def calculate_hydrophobicity(peptide):
    # Replace this with actual hydrophobicity calculation method.
    # Perhaps calculate the average hydrophobicity of amino acids in the peptide.
    # This is a dictionary that uses placeholder values. CAN CHANGE.
    hydrophobicity_scores = {
        'A': 0.17,
        'C': -0.24,
        'D': 1.23,
        # Add more amino acids and their hydrophobicity scores here.
    }
    # score = sum(hydrophobicity_scores.get(aa, 0) for aa in peptide)
    score = mean(hydrophobicity_scores.get(aa, 0) for aa in peptide)
    return score

# Calculate hydrophobicity scores for each peptide in the dataset
data['Hydrophobicity'] = data['Peptide'].apply(calculate_hydrophobicity)

# Display the modified dataset with added hydrophobicity scores
print(data.head())

# Split the dataset into two groups: binding peptides and non-binding peptides
binding_peptides = data[data['Bind'] == 1]['Hydrophobicity']
non_binding_peptides = data[data['Bind'] == 0]['Hydrophobicity']

# Perform a t-test to compare the hydrophobicity scores between the two groups
t_stat, p_value = stats.ttest_ind(binding_peptides, non_binding_peptides)

# Print the results
print("T-statistic:", t_stat)
print("P-value:", p_value)
# A T-statistic quantifies the size of the difference between groups
# A low p-value indicates a significant difference in hydrophobicity between binding and non-binding peptides.