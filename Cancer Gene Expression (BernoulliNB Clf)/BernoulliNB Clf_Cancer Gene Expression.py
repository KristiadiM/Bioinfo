"""Cancer Gene Expression (BernoulliNB).ipynb

# **Data Loading**
"""

import pandas as pd
import numpy as np

df = pd.read_csv('/content/Mod_Cancer_Dataset.zip')

"""# **Exploratory Data Analysis**"""

# Checking dataset size
print(f'Dataset shape: {df.shape}')

df.head()

# Checking columns
print(f'1st until 3rd columns: {df.columns[:3]}')
print(f'Last column: {df.columns[-1]}')

# Checking missing value
print(f'Is there any missing value? {df.isnull().values.any()}')
print('-'*30)

print(f'Null value: {df.isnull().sum()}')
print('-'*30)

# Show which column has missing value
missing_cols = df.columns[df.isnull().any()]
print("Columns with missing data:", list(missing_cols))

# Dataset type info
df.info()

print('-'*30)

df.dtypes

# Checking how many type of cancer
print(f'Cancer types: {df["Cancer Type"].unique()}')

# Visualize the cancer type distribution
import matplotlib.pyplot as plt

df['Cancer Type'].value_counts().plot(kind='bar')
plt.xlabel('Cancer Type')
plt.ylabel('Count')
plt.title('Cancer Type Distribution')
plt.show()

"""# **Data Preprocessing**"""

# Training set preparation
X = df.drop(columns=['Sample ID', 'Cancer Type'], axis=1)
y = df['Cancer Type']

X.shape

y.shape

"""**Features Engineering**"""

# Features Engineering with Label Encoding
from sklearn.preprocessing import LabelEncoder

label_encoder=LabelEncoder()
label_encoder.fit(y)
y_encoded=label_encoder.transform(y)
labels=label_encoder.classes_
classes=np.unique(y_encoded)

print("Original labels (Classes)):", labels)
print("Encoded labels:", classes)

"""**Data Splitting**"""

# Dataset Splitting
from sklearn.model_selection import train_test_split

X_train, X_test, y_train_encoded, y_test_encoded = train_test_split(X, y_encoded, test_size=0.2, random_state=42)

df.iloc[:,0:10].describe()

"""**Normalize Data**"""

# Data Normalization
from sklearn.preprocessing import MinMaxScaler

data_scaler = MinMaxScaler()
data_scaler.fit(X_train)
X_train = data_scaler.transform(X_train)
X_test = data_scaler.transform(X_test)

"""We have several features that have same value across all rows -> useless for classification

**Feature Selection**
"""

# Remove constant (or near-constant) features
from sklearn.feature_selection import VarianceThreshold

constant_filter = VarianceThreshold(threshold=0.0)
X_train = constant_filter.fit_transform(X_train)
X_test = constant_filter.transform(X_test)

# Check for NaNs
print("NaNs in X_train:", np.isnan(X_train).sum())

# Check for Infs
print("Infs in X_train:", np.isinf(X_train).sum())

# Select top K best features using ANOVA F-value
from sklearn.feature_selection import SelectKBest, f_classif

selector = SelectKBest(f_classif, k=500)
selector.fit(X_train, y_train_encoded)

# Apply transformation and keep only selected features
X_train_selected = selector.transform(X_train)
X_test_selected = selector.transform(X_test)
print("X_train_selected shape:", X_train_selected.shape)
print("X_test_selected shape:", X_test_selected.shape)

"""# **Model Selection and Training**

Using LazyPredict for Model Benchmarking
"""

!pip install lazypredict

import lazypredict
from lazypredict.Supervised import LazyClassifier

clf = LazyClassifier(verbose=0,ignore_warnings=True, custom_metric=None)
train_lazy,test_lazy = clf.fit(X_train_selected, X_test_selected, y_train_encoded, y_test_encoded)

print("LazyClassifier Train Results:")
train_lazy

print("LazyClassifier Test Results:")
test_lazy

"""We got BernoulliNB as the best algorithm (less time taken same result), we train our model with this algorithm

**Model Training**
"""

# Binarize the features
from sklearn.preprocessing import Binarizer

binarizer = Binarizer(threshold=0.0)
X_train_bin = binarizer.fit_transform(X_train_selected)
X_test_bin = binarizer.transform(X_test_selected)

"""**Model Evaluation**"""

# Model Training
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import classification_report

model = BernoulliNB()
model.fit(X_train_bin, y_train_encoded)

y_pred = model.predict(X_test_bin)
print(classification_report(y_test_encoded, y_pred, target_names=label_encoder.classes_))

"""# **Saving the Model**"""

import joblib

# Save the trained model to a file
joblib.dump(model, 'BernoulliNB_Cancer Predictor.pkl')

"""**Loading the Model**"""

# Load the saved model from the file
loaded_model = joblib.load('BernoulliNB_Cancer Predictor.pkl')
y_pred = loaded_model.predict(X_test_bin)

# Decoding the predicted labels back to their original class names
decoded_label = label_encoder.inverse_transform([y_pred[0]])
print(f"Predicted Cancer Type: {decoded_label[0]}")