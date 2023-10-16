import pandas as pd
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from scipy import stats
from statistics import mean, median
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
max_peptide_length = data['Peptide Length'].max()
print("Length of the longest peptide:", max_peptide_length)

amino_acid_list = 'ACDEFGHIKLMNPQRSTVWY'
hydrophobic = 'AVFPMIL' #1
polar = 'STYHCNQ' #2
charged = 'DEKR'#3
glycine = 'G' #4

average_position = []

# === FACTOR ENCODING ANALYSIS ===
# Define a function to encode amino acids based on a single factor
def factor_encoding(sequence):
    encoding = []
    positions = []
    position = 0
    counter = 0
    for aa in sequence:
        if (aa in polar):
            encoding.append('1')
            positions.append(position)
            counter+=1
        else:   
            encoding.append('0')
        position+=1
    while (len(encoding) < max_peptide_length):
        encoding.append('-2')
    if (len(positions) > 0):
        encoding.append(median(positions))
    else:
        encoding.append('-2')
    return encoding

# Add factor encoded column table for amino acid composition
amino_acid_features = data['Peptide'].apply(factor_encoding).tolist()
amino_acid_df = pd.DataFrame(amino_acid_features)
data = pd.concat([data, amino_acid_df], axis=1)

# Display the modified dataset with amino acid composition features
# print(data.head())

# Perform a logistic regression analysis
# Split the data into training and testing sets
X = data.drop(['Bind', 'Peptide', 'Peptide Length'], axis=1)  # Use one-hot encoded table as features
print(X)
y = data['Bind']  # Target variable

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a logistic regression classifier
# clf = DecisionTreeClassifier(random_state=42, max_depth=20)
clf = RandomForestClassifier(random_state=42, max_depth=20)
# clf = LogisticRegression(random_state=42)
clf.fit(X_train, y_train)

# Make predictions on the test set
y_pred = clf.predict(X_test)

# Evaluate the classifier
accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)
print("Classification Report:\n", classification_report(y_test, y_pred))
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

data.to_csv('output-polar.csv')