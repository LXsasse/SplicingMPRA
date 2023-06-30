import numpy as np
import sys, os
from sklearn import linear_model, ensemble
from scipy.stats import pearsonr
import matplotlib.pyplot as plt 

# check if number can be converted to a float
def isfloat(number):
    try:
        float(number)
    except:
        return False
    else:
        return True
    
# check if string can be integer or float
def numbertype(inbool):
    try:
        int(inbool)
    except:
        pass
    else:
        return int(inbool)
    try:
        float(inbool)
    except:
        pass
    else:
        return float(inbool)
    return inbool

# check if str is boolean, a list, or a number, otherwise return string back
import ast
def check(inbool):
    if inbool == 'True' or inbool == 'TRUE' or inbool == 'true':
        return True
    elif inbool == 'False' or inbool == 'FALSE' or inbool == 'false':
        return False
    elif inbool == 'None' or inbool == 'NONE' or inbool == 'none':
        return None
    elif "[" in inbool or "(" in inbool:
        return ast.literal_eval(inbool)
    else:
        inbool = numbertype(inbool)
    return inbool

# check if the string ftype is in one of the strings in files
def findkey(files, ftype):
    ptfi = []
    for i, fi in enumerate(files):
        if ftype in fi:
            ptfi.append(i)
    if len(ptfi) == 1:
        return files[ptfi[0]]
    else:
        return None

# Read in values for X and Y, sort and align them, and return them with names of data points, output tracks 
def readin(inputfile, outputfile, delimiter = ' ', return_header = True, skip_header = 0):
    # inputfile is a npz file
    Xin = np.load(inputfile, allow_pickle = True)
    xfiles = Xin.files
    fname =findkey(xfiles,'features') 
    # In earlier version, the feature matrix was combined with the feature names
    X = Xin[fname]
    if len(X) == 2:
        X, inputfeatures = X
    else:
        if 'featuretypes' in Xin.files:
            inputfeatures = Xin['featuretypes']
        else:
            inputfeatures = np.arange(np.shape(X)[-1])
    dname =findkey(xfiles,'names')
    inputnames = Xin[dname]
    
    # file with target values can be npz or txt file
    if os.path.isfile(outputfile):
        if os.path.splitext(outputfile)[1] == '.npz':
            Yin = np.load(outputfile, allow_pickle = True)
            Y, outputnames = Yin['counts'], Yin['names'] # Y should of shape (nexamples, nclasses, l_seq/n_resolution)
        else:
            Yin = np.genfromtxt(outputfile, dtype = str, delimiter = delimiter, skip_header = skip_header)
            Y, outputnames = Yin[:, 1:].astype(float), Yin[:,0]
        hasoutput = True
    else:
        print(outputfile, 'not a file')
        hasoutput = False
        Y, outputnames = None, None
    
    # sort X and Y so that data points align
    sortx = np.argsort(inputnames)
    if hasoutput:
        sortx = sortx[np.isin(np.sort(inputnames), outputnames)]
        sorty = np.argsort(outputnames)[np.isin(np.sort(outputnames), inputnames)]
        X, inputnames = X[sortx], inputnames[sortx]
        outputnames = outputnames[sorty]
        Y = Y[sorty]
        if not np.array_equal(inputnames, outputnames):
            print('Sorting was not successful, input and output do not align')
            sys.exit()
    
    # remove non-existant features
    mask = np.sum(np.abs(X), axis = 0)
    X, inputfeatures = X[:,mask], inputfeatures[mask]
    
    # if several target tracks are in target file, also provide name in header if given in file
    if return_header and hasoutput:
        if os.path.splitext(outputfile)[1] == '.npz':
            if findkey(Yin.files, 'columns') is not None:
                header = Yin[findkey(Yin.files, 'columns')]
            else:
                header = ['C'+str(i) for i in range(np.shape(Y)[1])]
        else:
            header = open(outputfile, 'r').readline()
            if '#' in header:
                header = header.strip('#').strip().split(delimiter)
            else:
                header = ['C'+str(i) for i in range(np.shape(Y)[1])]
        header = np.array(header)
    else:
        header  = None
    
    print('Input shapes X', np.shape(X))
    print('Output shapes Y', np.shape(Y))
        
    return X, Y, inputnames, inputfeatures, header

## Create random train, test and validation sets using permutations
# If samples can be split into different classes, each set will have same fraction of each class
def create_sets(n_samples, folds, fold, seed = 1010, Yclass = None, genenames = None):
    if isinstance(folds, int):
        if Yclass is None:
            Yclass = np.ones(n_samples)
        np.random.seed(seed)
        permut = np.random.permutation(n_samples)
        classes = np.unique(Yclass)
        testset, valset, trainset = [],[],[]
        valfold = fold -1
        if valfold == -1:
            valfold = folds -1
        for cla in classes:
            inclass = Yclass[permut] == cla
            totinclass = np.sum(inclass)
            modulos = totinclass%folds
            addedstart = fold * int(fold<modulos)
            startfold = int((fold/folds)*np.sum(Yclass== cla)) +addedstart
            endfold = int(((fold+1)/folds)*np.sum(Yclass== cla))+ addedstart + int(fold+1 < modulos)
            valaddedstart = valfold * int(valfold<modulos)
            valstartfold = int((valfold/folds)*np.sum(Yclass== cla)) +valaddedstart
            valendfold = int(((valfold+1)/folds)*np.sum(Yclass== cla))+ valaddedstart + int(valfold+1 < modulos)
            testset.append(permut[inclass][startfold:endfold])
            valset.append(permut[inclass][valstartfold:valendfold])
        testset, valset = np.concatenate(testset), np.concatenate(valset)
        trainset = np.delete(np.arange(n_samples, dtype = int), np.append(testset, valset))
    elif os.path.isfile(folds):
        lines = open(folds, 'r').readlines()
        sets = []
        for l, line in enumerate(lines):
            if line[0] != '#':
                line = line.strip().split()
                sets.append(np.where(np.isin(genenames,line))[0])
        testset = sets[(fold+1)%len(sets)]
        valset = sets[fold]
        trainset = np.delete(np.arange(len(genenames), dtype = int),np.append(testset,valset))
    return trainset, testset, valset

# normalize X
def manipulate_input(X, features, sysargv, outname):
    if '--select_features' in sysargv:
        selected_feat = np.genfromtxt(sysargv[sysargv.index('--select_features')+1], dtype = str)
        featmask = np.isin(features, selected_feat)
        features, X = np.array(features)[featmask], X[:, featmask]
        print('X reduced to', np.shape(X))
        outname+= '_featsel'+str(len(X[0]))
        
    if '--centerfeature' in sysargv:
        outname += '-cenfeat'
        X = X - np.mean(X, axis = 0)
    
    elif '--centerdata' in sysargv:
        outname += '-cendata'
        X = X - np.mean(X, axis = 1)[:, None]

    if '--norm2feature' in sysargv:
        outname += '-n2feat'
        norm = np.sqrt(np.sum(X*X, axis = 0))
        norm[norm == 0] = 1.
        X = X/norm
        
    elif '--norm1feature' in sysargv:
        outname += '-n1feat'
        norm =np.sum(np.absolute(X), axis = 0)
        norm[norm == 0] = 1.
        X = X/norm
    
    if '--norm2data' in sysargv:
        outname += '-n2data'
        norm =np.sqrt(np.sum(X*X, axis = 1))[:, None] 
        X = X/norm
        
    elif '--norm1data' in sysargv:
        outname += '-n1data'
        X = X/np.sum(np.absolute(X), axis = 1)[:,None]
    
    if '--standardfeature' in sysargv:
        outname += '-standfeat'
        norm = np.sqrt(np.mean(X*X, axis = 0))
        norm[norm == 0] = 1.
        X = (X-np.mean(X,axis = 0))/norm
        
    if '--standarddata' in sysargv:
        outname += '-standdata'
        norm =np.sqrt(np.mean(X*X, axis = 1))[:, None] 
        X = (X-np.mean(X, axis = 1)[:,None])/norm
    
    return X, features, outname

        
def mse(x,y):
    return (1./np.prod(np.shape(x)))* np.sqrt(np.sum((x-y)**2))

# Save performance of method
def performance(Y_pred, Ytest, outname, sysargv, compare_random = True):
    mses = mse(Ytest, Ypred)
    corr, pcorr = pearsonr(Ytest, Ypred)
    obj = open(outname+'_performance.txt', 'w')
    obj.write('MSE\t'+str(mses)+'\n')
    obj.write('Pearson R\t'+str(corr)+'\n')
    obj.write('Pearson R pvalue\t'+str(pcorr)+'\n')
    print('MSE\t'+str(mses)+'\nPearson R\t'+str(corr)+'\n'+'Pearson R pvalue\t'+str(pcorr)+'\n')
    if compare_random:
        randmse = mse(Ytest, np.random.permutation(Ypred))
        randcorr = pearsonr(Ytest, np.random.permutation(Ypred))[0]
        print('Random MSE\t'+str(randmse)+'\nRandom Pearson R\t'+str(randcorr)+'\n')

# Make a scatter plot
def plot_scatter(Ytest, Ypred, titles = None, xlabel = None, ylabel = None, outname = None, include_lr = True, include_mainvar = True):
    n = len(Ytest[0])
    if n > 100:
        print('Number of examples is too large', n)
        return
    x_col = int(np.sqrt(n))
    y_row = int(n/x_col) + int(n%x_col!= 0)
    fig = plt.figure(figsize = (x_col*3.5,y_row*3.5), dpi = 100)
    for e in range(n):
        ax = fig.add_subplot(y_row, x_col, e+1)
        pcorr = pearsonr(Ytest[:,e], Ypred[:,e])[0]
        if titles is not None:
            ax.set_title(titles[e]+' R='+str(np.around(pcorr,2)), fontsize = 6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.scatter(Ytest[:,e], Ypred[:,e], c = 'slategrey', alpha = 0.7, s = 6)
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        limx, limy = ax.get_xlim(), ax.get_ylim()
        lim = [max(limx[0], limy[0]), min(limx[1], limy[1])]
        ax.plot(lim, lim, color = 'dimgrey', ls = '-')
        if include_lr:
            lr = linear_model.LinearRegression().fit(Ytest[:, [e]], Ypred[:,e])
            ax.plot(np.array(limx), lr.predict(np.array(limx).reshape(-1,1)), color = 'r')
        if include_mainvar:
            centerx, centery = np.mean(Ytest[:,e]), np.mean(Ypred[:,e])
            maindir, u, v = np.linalg.svd(np.array([Ytest[:,e]-centerx, Ypred[:,e]-centery]), full_matrices=False)
            maindir = maindir[:,0]
            slope = maindir[1]/maindir[0]
            bias = centery-slope*centerx
            ax.plot(np.array(limx), np.array(limx)*slope + bias, color = 'r')
        
    if xlabel is not None:
        fig.text(0.5, 0.05-0.25/y_row, xlabel, ha='center')
    if ylabel is not None:
        fig.text(0.05-0.2/x_col, 0.5, ylabel, va='center', rotation='vertical')
    if outname is not None:
        print('SAVED as', outname)
        fig.savefig(outname, dpi = 200, bbox_inches = 'tight')
    else:
        fig.tight_layout()
        plt.show()

# define the predictor to be used        
class feature_predictor():
    def __init__(self, model, **params):
        self.model = model
        self.params = params
        
        if model == 'RandomForest' or model == 'RF' or model == 'randomforest':
            self.lr = ensemble.RandomForestRegressor(**params) #n_estimators=100, *, criterion='squared_error', max_depth=None, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_features=1.0, max_leaf_nodes=None, min_impurity_decrease=0.0, bootstrap=True, oob_score=False, n_jobs=None, random_state=None, verbose=0, warm_start=False, ccp_alpha=0.0, max_samples=None)
        
        if model == 'LinearRegression' or model == 'OLS' or model == 'LR':
            self.lr = linear_model.LinearRegression(**params) #, fit_intercept=True, copy_X=True, n_jobs=None, positive=False)
        
        if model == 'ElasticNet' or model == 'elastic' or model == 'EN' or model == 'elasticnet':
            self.lr = linear_model.ElasticNet(**params) #alpha=1.0, *, l1_ratio=0.5, fit_intercept=True, precompute=False, max_iter=1000, copy_X=True, tol=0.0001, warm_start=False, positive=False, random_state=None, selection='cyclic')
        
    def fit(self, x, y):
        self.lr=self.lr.fit(x,y)
        self.x = x
        self.y = y
        
    def predict(self, x):
        pred = self.lr.predict(x)
        return pred
    
    def feature_importance(self):
        # assigns t-test statistic to LR coefficients, same as statsmodels.api.OLS
        if self.model == 'LinearRegression' or self.model == 'OLS' or self.model == 'LR' or self.model == 'ElasticNet' or self.model == 'elastic' or self.model == 'EN' or self.model == 'elasticnet':
            if self.lr.fit_intercept:
                params = np.append(self.lr.intercept_, self.lr.coef_)
                newX = np.append(np.ones((len(self.x),1)), self.x, axis=1)
            else:
                params = self.lr.coef_
                newX = self.x
            predictions = self.lr.predict(self.x)
            MSE = (sum((self.y-predictions)**2))/(len(newX)-len(newX[0]))
            var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
            sd_b = np.sqrt(var_b)
            ts_b = params/sd_b
            #p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b]
            
            '''
            import statsmodels.api as sm
            X2 = sm.add_constant(self.x)
            est = sm.OLS(self.y, X2)
            est2 = est.fit()
            print(est2.summary())
            print(ts_b)
            print(self.lr.coef_)
            sys.exit()
            '''
            
            if self.lr.fit_intercept:
                ts_b = ts_b[1:]
            rank = np.argsort(np.argsort(-np.abs(ts_b))) + 1
            return rank
        
        elif self.model == 'RandomForest' or self.model == 'RF' or self.model == 'randomforest':
            ftimp = self.lr.feature_importances_
            return np.argsort(np.argsort(-ftimp)) + 1
    
    def coef_(self):
        if self.model == 'RandomForest' or self.model == 'RF' or self.model == 'randomforest':
            return self.lr.feature_importances_
        return self.lr.coef_
    
    def intercept_(self):
        if self.model == 'RandomForest' or self.model == 'RF' or self.model == 'randomforest':
            return 0
        return self.lr.intercept_
  
# Save the feature statistics
def feature_stats(model, featurenames, outname):
    featrank = model.feature_importance()
    feat_coef = model.coef_()
    sort = np.argsort(featrank)
    for c, co in enumerate(sort[:20]):
        print(featurenames[co], feat_coef[co], featrank[co])
    np.savetxt(outname+'_features.txt', np.array([featurenames, feat_coef, featrank], dtype = str).T, fmt = '%s')
  
# Combine filenames to a new output file name, removing text that is redundant in both filenames    
def create_outname(name1, name2, lword = 'on'):
    name1s = os.path.split(name1)[1].replace('.dat','').replace('.hmot','').replace('.txt','').replace('.npz','').replace('.list','').replace('.csv','').replace('.tsv','').replace('.tab','').replace('-', "_").replace('.fasta','').replace('.fa','').split('_')
    name2s = name2.replace('-', "_").split('_')
    diffmask = np.ones(len(name1s)) == 1
    for n, na1 in enumerate(name1s):
        for m, na2 in enumerate(name2s):
            if na1 in na2:
                diffmask[n] = False
    diff = np.array(name1s)[diffmask]
    outname = os.path.split(name2)[1].replace('.dat','').replace('.hmot','').replace('.txt','').replace('.npz','').replace('.list','').replace('.csv','').replace('.tsv','').replace('.tab','').replace('.fasta','').replace('.fa','')
    if len(diff) > 0:
        outname+=lword+'_'.join(diff)
    return outname








if __name__ == '__main__':
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
    
    outname = create_outname(inputfile, outputfile)
    
    delimiter = ' ' 
    if '--delimiter' in sys.argv:
        delimiter = sys.argv[sys.argv.index('--delimiter')+1]
    
    skip_header = 0
    if '--skip_header' in sys.argv:
        skip_header = int(sys.argv[sys.argv.index('--skip_header')+1])
    
    # read in data files
    X, Y, datanames, featurenames, header = readin(inputfile, outputfile, delimiter = delimiter, return_header = True, skip_header = skip_header)
    # Normalize input features
    if np.shape(Y)[1] == 1:
        Y = Y.reshape(np.shape(Y)[0])
    X, features, outname = manipulate_input(X, featurenames, sys.argv, outname)
    
    folds, fold = 10,0
    if '--crossvalidation' in sys.argv:
        folds, fold = numbertype(sys.argv[sys.argv.index('--crossvalidation') + 1]), numbertype(sys.argv[sys.argv.index('--crossvalidation') + 2])
        outname += 'fold'+str(fold)
    # create training and test sets
    trainset, testset, valset= create_sets(len(X), folds, fold, seed = 1010, Yclass = None, genenames = datanames)
    
    if '--combine_train_and_val' in sys.argv:
        trainset = np.append(trainset, valset)
    else:
        testset = valset
        
    # define model specifics
    if '--model_params' in sys.argv:
        modeltype = sys.argv[sys.argv.index('--model_params')+1]
        outname += modeltype
        params = {}
        if len(sys.argv) > sys.argv.index('--model_params')+2:
            if '--' not in sys.argv[sys.argv.index('--model_params')+2]:
                if '+' in sys.argv[sys.argv.index('--model_params')+2]:
                    parameters = sys.argv[sys.argv.index('--model_params')+2].split('+')
                else:
                    parameters = [sys.argv[sys.argv.index('--model_params')+2]]
                for p in parameters:
                    if ':' in p and '=' in p:
                        p = p.split('=',1)
                    elif ':' in p:
                        p = p.split(':',1)
                    elif '=' in p:
                        p = p.split('=',1)
                    params[p[0]] = check(p[1])
                    outname += p[0][:2]+p[0][max(2,len(p[0])-2):]+str(p[1])
        model = feature_predictor(modeltype, **params)
    
    model.fit(X[trainset],Y[trainset])
    Ypred = model.predict(X[testset])
    
    performance(Ypred, Y[testset], outname, sys.argv) 
    
    feature_stats(model, featurenames, outname)
 
    plot_scatter(Y[testset].reshape(-1,1), Ypred.reshape(-1,1), titles = [''], xlabel = 'Measured', ylabel = 'Predicted '+modeltype, outname = outname + '.jpg', include_lr = True, include_mainvar = False)
    
    
    



