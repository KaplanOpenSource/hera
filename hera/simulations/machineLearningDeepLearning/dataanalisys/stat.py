import numpy as np
import sklearn.metrics

def stat(a,b, kind='rmse', nantozero=True, verbose=None):
    # calculating the score of comparing two 1D vectors
    # kind can be: rmse, mae, bias, r2, r, nmse, fb
    # nantozero, if True, convert nan values to 0, 
    #            if False, comparing the lists without the nans, if the length is not the same, use only the first ones.

    # return np.sqrt(np.mean((predictions-targets)**2))
    # np.corrcoef(db[i][sel], db[-1][sel])[1, 0]
    try:
        if nantozero:
            a[np.isnan(a)]=0.0
            b[np.isnan(b)]=0.0
        else:
            # not implemented yet, but if we remove nans, we should remove the exact cells in the other list
            a=a[~np.isnan(a)]
            b=b[~np.isnan(b)]
            minlen=min(len(a),len(b))
            a=a[:minlen]
            b=b[:minlen]
        if (len(a)<1) or (len(b)<1):
            c=-999
        else:
            if kind=='rmse':
                c = np.sqrt(np.mean((a-b)**2))
            elif kind=='mae':
                c = np.mean(np.abs(a-b))
            elif kind=='bias':
                c = np.mean(a-b)
            elif kind=='r2':
                # coefficient of determination
                c = sklearn.metrics.r2_score(a.ravel(), b.ravel())
            elif kind=='r':
                # Pearson correlation coefficient 
                c = np.corrcoef(a, b)[1, 0]
        #        if len(c.shape)>1:
        #            c=c[1,0]
            elif kind=='nmse':
                # normalized mean square error, chang and hanna 2004
                c = np.sqrt(np.mean((a-b)**2))/(np.mean(a)*np.mean(b))
            elif kind=='fb':
                # Fractional bias, chang and hanna 2004
                c = (np.mean(a)-np.mean(b))/(0.5*(np.mean(a)+np.mean(b)))
            else:
                print('unknown kind ',kind)
                c=-998
    except:
        c=-997
        print('stat except', a, b, len(a), len(b), kind)

    return round(c,4)

if __name__ == '__main__':
    a=[1,2,3,4]
    b=[2,3,4,7]
    print(stat(np.asarray(a),np.asarray(b),kind='rmse'))
