# Linear Support Vector Classifier with L1 regularization.
# Requires inputs X and y already created.  Vary C parameter.


C_min = sklearn.svm.l1_min_c(X, y, fit_intercept=True, intercept_scaling=1.0)
SVC = sklearn.svm.LinearSVC(C=100*C_min, penalty='l1', dual=False)


# K-fold crossvalidation and ROC curve                                                                                                                                                                             
N_CV_repeats = 10
mean_tpr = 0.0                                                                                                                                                                                                      
mean_fpr = np.linspace(0, 1, 100) 
tpr_vals = []
auc_vals = []

for N in range(N_CV_repeats):

    cv = StratifiedKFold(y, 3, shuffle=True)
    for train_index, test_index in cv:
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index] 
        y_score = SVC.fit(X_train, y_train).decision_function(X_test)

        # Compute ROC curve and ROC area for each class
        fpr, tpr, thresholds = roc_curve(y_test, y_score, pos_label='1')
        roc_auc = auc(fpr, tpr)
        auc_vals.append(roc_auc)
        counter += 1
        tpr_interpolated = scipy.interp(mean_fpr, fpr, tpr)
        tpr_vals.append(tpr_interpolated)
        mean_tpr += tpr_interpolated
    
mean_tpr /= len(cv)*(N_CV_repeats)
mean_auc = np.mean(auc_vals)
print "" 
print "Random Forest AUC = " + str(np.mean(auc_vals)) + ' +/- ' + str(np.std(auc_vals))
print ""

# Calculate Mann-Whitney U statistic and p-value
n1 = np.sum(y == '1')
n2 = np.sum(y == '0')
U = mean_auc * n1 * n2
mu = n1*n2/2
sigma = np.sqrt(n1*n2*(n1+n2+1)/12)
Z_score_MWU = (U - mu)/sigma
pval = 1-scipy.stats.norm.cdf(Z_score_MWU)
print "Z-score for AUC (Mann-Whitney U test): " + str(Z_score_MWU) 
print "p-value (assuming Gaussian) = " + str(pval)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
sns.set_style('white')
plt.plot(mean_fpr, mean_tpr, color='black', lw=3) 
plt.plot(mean_fpr, mean_fpr, '-r')
plt.title('AUC = ' + str(mean_auc)[:5], fontsize=20)
std_tpr = np.std(tpr_vals, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')    
sns.despine()
plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel('False positive rate', fontsize = 25)
plt.ylabel('True positive rate', fontsize = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)


