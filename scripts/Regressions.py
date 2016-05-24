"""

OVERVIEW: 

Module for automatically performing regressions on provided data.

"""

import numpy as np
import pandas
import sklearn.linear_model
import sklearn.ensemble
import matplotlib
import matplotlib.pyplot as plt

def compute_linear_regressions(features, metadata):
    # Computes linear regressions of the feature vectors against each metadata type.
    # Takes as input a feature object (from Features.py) module, and a metadata dictionary
    # which has the form 'metadata_dict[sampleID][metadata_type]', where sample IDs are the same
    # in both feature object and metadata dictionary.
    regr = sklearn.linear_model.LinearRegression()

    # Initialize
    nFeatures = len(features.feature_vectors[features.feature_vectors.keys()[0]])
    nSamples = len(features.feature_vectors.keys())
    sampleIDs = features.feature_vectors.keys()
    nMetadata = len(metadata[metadata.keys()[0]])
    metadataTypes = [key for key in metadata[metadata.keys()[0]].keys()]

    thetas = np.zeros((nMetadata, nFeatures))
    intercepts = np.zeros((nMetadata,1))
    mean_square_errors = np.zeros((nMetadata,1))
    Rsquared_vals = np.zeros((nMetadata,1))    

    # Prepare input data matrix (A, for Ax = b)
    A = np.zeros((nSamples, nFeatures)) 
    for i in range(nSamples):
        A[i,:] = features.feature_vectors[sampleIDs[i]]

    # Loop through each metadata type and perform a linear regression
    for i in range(len(metadataTypes)):
        b = np.zeros((nSamples, 1))
        for (index, sampleID) in enumerate(sampleIDs):
            b[index,0] = metadata[sampleID][i] 
            regr.fit(A, b)
            thetas[i,:] = regr.coef_
            intercepts[i,0] = regr.intercept_
            Rsquared_vals[i] = regr.score(A,b)
    
            # Record the mean square error, as well as this error normalised by mean abundance
            mean_square_errors[i,0] = np.mean((regr.predict(A) - b) ** 2)

    return [thetas, intercepts, Rsquared_vals, mean_square_errors] 
    
def compute_linear_regressions_lasso(features, metadata, iterations=1000):
    # Computes linear regressions of the feature vectors against each metadata type.
    # Takes as input a feature object (from Features.py) module, and a metadata dictionary
    # which has the form 'metadata_dict[sampleID][metadata_type]', where sample IDs are the same
    # in both feature object and metadata dictionary.
    clf = sklearn.linear_model.Lasso(alpha=0.1, max_iter=iterations)
    
    # Initialize
    nFeatures = len(features.feature_vectors[features.feature_vectors.keys()[0]])
    nSamples = len(features.feature_vectors.keys())
    sampleIDs = features.feature_vectors.keys()
    nMetadata = len(metadata[metadata.keys()[0]])
    metadataTypes = [key for key in metadata[metadata.keys()[0]].keys()]

    thetas = np.zeros((nMetadata, nFeatures))
    intercepts = np.zeros((nMetadata,1))
    mean_square_errors = np.zeros((nMetadata,1))
    Rsquared_vals = np.zeros((nMetadata,1))    

    # Prepare input data matrix (A, for Ax = b)
    A = np.zeros((nSamples, nFeatures)) 
    for i in range(nSamples):
        A[i,:] = features.feature_vectors[sampleIDs[i]]

    # Loop through each metadata type and perform a linear regression
    for i in range(len(metadataTypes)):
        b = np.zeros((nSamples, 1))
        for (index, sampleID) in enumerate(sampleIDs):
            b[index,0] = metadata[sampleID][i] 
            clf.fit(A, b)
            thetas[i,:] = clf.coef_
            intercepts[i,0] = clf.intercept_
            Rsquared_vals[i] = clf.score(A,b)
            
    
            # Record the mean square error, as well as this error normalised by mean abundance
            mean_square_errors[i,0] = np.mean((clf.predict(A) - b) ** 2)

    return mean_square_errors

def compute_random_forest_regressions(features, metadata, figure_name):
    # Takes as input a feature object and a metadata dictionary (as in the regressions function above).
    # Trains an RF regressor against each metadtata type.

    # Build RF regressor
    RFregr = sklearn.ensemble.RandomForestRegressor(n_estimators=1000)

    # Initialize
    nFeatures = len(features.feature_vectors[features.feature_vectors.keys()[0]])
    nSamples = len(features.feature_vectors.keys())
    sampleIDs = features.feature_vectors.keys()
    nMetadata = len(metadata[metadata.keys()[0]])
    metadataTypes = [key for key in metadata[metadata.keys()[0]].keys()]

    # Prepare input data matrix (A, for Ax = b)
    X = np.zeros((nSamples, nFeatures)) 
    for i in range(nSamples):
        X[i,:] = features.feature_vectors[sampleIDs[i]]
    
    # Loop through each metadata type and train an RF regressor
    f, axarr = plt.subplots(0,nMetadata, sharex=True)
    axes_counter = 0
    for i in range(len(metadataTypes)):
        y = np.zeros((nSamples, 1))
        for (index, sampleID) in enumerate(sampleIDs):        
            y[index,0] = metadata[sampleID] 
        
        # Train RF regressor
        cv = StratifiedKFold(y, 7)
        counter = 1
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        for train_index, test_index in cv:
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            probabs = diet_RF.fit(X_train, y_train).predict_proba(X_test)
            # Compute ROC curve and area the curve
            fpr, tpr, thresholds = roc_curve(y_test, probabs[:, 1], pos_label=cluster_labels[key])
            counter += 1        
            mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr /= len(cv)
        roc_auc = auc(mean_fpr, mean_tpr)
        axarr[0,axes_counter].plot(mean_fpr, mean_tpr, lw=1)
        axarr[0,axes_counter].set_title(key)
        axarr[0,axes_counter].axis([0, 1, 0, 1])
        axes_counter += 1
        print "AUC for " + key + ": " + str(roc_auc)
    savefig(figure_name)

    
