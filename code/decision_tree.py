from sklearn.tree import DecisionTreeClassifier
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
from dtreeviz.trees import dtreeviz
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score
data = pd.read_csv('~/Documents/MBP/TaboriLab/codes/mutation_project/gene_mutation_project/vcf_validation/output/orientation_counts/orientation_counts_normality_test.csv')
WRITE = True
#replace good with 1 and bad with 0
data['tmb_guess'] = data['tmb_guess'].replace(['good','bad'],[1,0])

#train test split
X_train, X_test, y_train, y_test = train_test_split(data.iloc[:,:-1], data.iloc[:,-1],test_size=0.2,stratify=data.iloc[:,-1],random_state=42)


dt = DecisionTreeClassifier(random_state=0,max_depth=4,min_samples_leaf=3,criterion='gini',splitter='random').fit(X_train.iloc[:,1:], y_train)
print(np.average(cross_val_score(dt, X_train.iloc[:,1:], y_train, cv=10)))
print(accuracy_score(y_test, dt.predict(X_test.iloc[:,1:])))

viz = dtreeviz(dt, 
    np.array(X_train.iloc[:,1:]), 
    np.array(y_train),
    target_name='tmb_guess',
    feature_names=data.columns[1:-1], 
    class_names=['bad','good']
    ) 

viz.view()