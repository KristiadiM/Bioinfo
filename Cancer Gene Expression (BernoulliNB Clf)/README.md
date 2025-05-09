# Cancer Gene Expression Classification using Bernoulli Naive Bayes
This project implements a machine learning model to classify cancer types based on gene expression data. 
The dataset contains gene expression levels for various samples, and the goal is to predict the type of cancer using these features.

## Requirements:
pip install pandas numpy matplotlib scikit-learn lazypredict joblib

## Dataset Overview
- Dataset Source : Mod_Cancer_Dataset.zip (modified from https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq, merged the data and labels)
- Features : Gene expression levels (continuous numerical values)
- Target Variable : Cancer Type (multi-class classification)
- Sample ID : Unique identifier for each sample
- Data Size : (n_samples, n_features) — dynamically shown during execution

## Exploratory Data Analysis (EDA)
- Missing Values : Checked for nulls; none found.
- Cancer Types : Multiple cancer types identified in the dataset (e.g., Breast, Lung, Prostate).
- Distribution Visualization : Bar chart showing class distribution for cancer types.
  
## Data Preprocessing
Steps:
- Dropped Irrelevant Columns :
- Removed Sample ID as it doesn't contribute to prediction.
- Label Encoding :
    Encoded the target labels (Cancer Type) into numeric classes.
- Train-Test Split :
    80% training and 20% testing split with fixed random seed (random_state=42).
- Normalization :
    Used MinMaxScaler to scale feature values between 0 and 1.
- Feature Selection :
    Removed constant or near-constant features using VarianceThreshold.
- Selected top 500 features using ANOVA F-value (SelectKBest).

## Model Selection with LazyPredict
Used LazyClassifier from the lazypredict library to benchmark multiple classifiers. Results showed that BernoulliNB achieved competitive accuracy with minimal computational cost.

## Model Training: Bernoulli Naive Bayes
Implementation Details:
- Binarization : Applied binarizer thresholding (threshold=0.0) to convert continuous features into binary inputs suitable for BernoulliNB.
- Model : Trained using BernoulliNB() classifier.
Evaluation Metrics :
- Precision
- Recall
- F1-score
- Accuracy
- Confusion matrix (implied via classification report)
  
## Evaluation Report
The final model was evaluated using a classification report, which includes precision, recall, and F1-scores per class, along with overall accuracy.

## Saving and Loading the Model
- Saved Model : Using joblib.dump(), the trained model is saved as BernoulliNB_Cancer Predictor.pkl.
- Loading Model : The saved model can be reloaded using joblib.load() for future predictions.
