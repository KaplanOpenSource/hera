import numpy as np
import datetime
import matplotlib.pyplot as plt
import pickle

from sklearn.linear_model import LinearRegression, Lasso, ElasticNet, LogisticRegression, SGDRegressor, Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor

from sklearn.svm import LinearSVR
from sklearn import svm
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.decomposition import PCA

from sklearn import preprocessing
import sklearn.metrics

#from sklearn.cross_decomposition import CCA

        
import multiprocessing
import time

#model = SVR(cache_size=7000)

class ml(object):

    """
    A class that wrap the machine learning options

    """
    def __init__(self):
        """
        Choosing which ml kernel to use
        Currently we     have only one scaler
        """
        self.models = []
#        self.models.append(['SGD', SGDRegressor(max_iter=1000)])
  #        self.models.append(['Lasso', Lasso(alpha=0.1)])
  #        self.models.append(['SGDRegressor', SGDRegressor()])
#        self.models.append(['ElasticNet', ElasticNet(random_state=0)])
#        self.models.append(['GradientBoostingRegressor', GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=0, loss='ls')])
#        self.models.append(['LR', LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')])
#        self.models.append(['KNN2', KNeighborsRegressor(n_neighbors=2, weights = 'distance')])
#        self.models.append(['KNN4', KNeighborsRegressor(n_neighbors=4)])
        self.models.append(['KNN6', KNeighborsRegressor(n_neighbors=6, weights = 'distance')])
#        self.models.append(['KNN8', KNeighborsRegressor(n_neighbors=8)])
#        self.models.append(['KNN10', KNeighborsRegressor(n_neighbors=10)])
#        self.models.append(['KNN12', KNeighborsRegressor(n_neighbors=12)])
#        self.models.append(['SVM', svm.SVR(C=0.001, tol=1e-01)])
#        self.models.append(['LinratSVR', LinearSVR(C=1.0, dual=True, epsilon=0.0, fit_intercept=True, intercept_scaling=1.0, loss='epsilon_insensitive', max_iter=1000, random_state=0, tol=1e-05, verbose=0)])
         #self.models.append(['Ridge', Ridge(alpha=1.0)])
#        self.models.append(['LR', LinearRegression()])
#        self.models.append(['ANN', MLPRegressor(random_state=1, max_iter=500)])
        # self.models.append(['RF', RandomForestRegressor(n_estimators = 100, random_state = 0)])
        # self.models.append(['XGB', MLPRegressor(random_state=1, max_iter=500)])
        print('normalized data')

    def save(self, filename):
        print('mlsave',filename)
        with open(filename+'model.pkl', 'wb') as output:
            pickle.dump(self.model, output, pickle.HIGHEST_PROTOCOL) 
        print ('mid ml save')
        with open(filename+'scaler.pkl', 'wb') as output:
            pickle.dump(self.scaler, output, pickle.HIGHEST_PROTOCOL)             
        print('saved ml')
    
    def load(self, filename):
        with open(filename+'model.pkl', 'rb') as input:
            self.model = pickle.load(input)
        print('model;', self.model)
        with open(filename+'scaler.pkl', 'rb') as input:
            self.scaler = pickle.load(input)
        print('scaler;', self.scaler)
        
    def features_coeff(self, features, labels, featurestest = None, labeltest = None):
#       https://stats.stackexchange.com/questions/311488/summing-feature-importance-in-scikit-learn-for-a-set-of-features
        
        msklabel = alllabels == 0.0 #remove all for 3 lines
        features0 = allfeatures[~msklabel]
        labels0 = alllabels[~msklabel]

        if featurestest is None:
            print ('random sample')
            self.msk = np.random.rand(len(labels0)) < 0.9
            trainx = features0[self.msk]
            trainy = labels0[self.msk]
            testx = features0[~self.msk]
            testy = labels0[~self.msk]
        else:
            print ('chosen sample')
            msktest = labels == 0.0
            features1 = features[~msktest]
            labels1 = labels[~msktest]
            trainx = features0
            trainy = labels0
            testx = features1
            testy = labels1

            
        scaler = preprocessing.MinMaxScaler()
        scaler.fit(trainx)
        strainx = scaler.transform(trainx)
#        print('strainx', strainx, datetime.datetime.now())
        from sklearn.ensemble import RandomForestClassifier
#                rf = RandomForestClassifier(random_state=42).fit(features0, labelsux0)
        from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
        reg = LassoCV()
        reg.fit(strainx,trainy)
#                print(rf.feature_importances_)
        print(reg.coef_)
                
            
    def fit(self, features, labels, show='', featurestest = None, labeltest = None):

        if 4==5:
            models = [] 
            models.append(['KNN6', KNeighborsRegressor(n_neighbors=6, weights = 'distance')])
#            alllabels=np.load('nbp.npy').ravel()

            alllabelsu=np.load('nbu.npy').ravel()
            alllabelsv=np.load('nbv.npy').ravel()
            alllabelsw=np.load('nbw.npy').ravel()
            alllabels = (alllabelsu**2+alllabelsv**2+alllabelsw**2)**0.5
#            alllabels = np.arctan2(alllabelsv,alllabelsu)*180/np.pi
            jpjp=500
            alllabelsu=alllabelsu[::jpjp]
            alllabelsv=alllabelsv[::jpjp]
            alllabelsw=alllabelsw[::jpjp]
            # del alllabelsu, alllabelsv, alllabelsw
            allfeatures=np.load('nbaf.npy')
            grid_y=np.load('nbgrid_y.npy').ravel()
            alllabels=alllabels[::jpjp]
            allfeatures=allfeatures[::jpjp]
            grid_y=grid_y[::jpjp]
            plt.figure()
            plt.plot(alllabels,'.')
            
            miny = grid_y.ravel().min()
            maxy = grid_y.ravel().max()
            part1y = miny + (maxy-miny) * 0.25# 25
            part2y = miny + (maxy-miny) * 0.5
            part1y1 = miny + (maxy-miny) * 0.5# 25
            part2y1 = miny + (maxy-miny) * 0.75
            features = allfeatures[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            featurestest = allfeatures[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]
            labels = alllabels[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsu = alllabelsu[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsv = alllabelsv[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelslimit = alllabelsu[(grid_y.ravel()>part1y) & (grid_y.ravel()<part2y)]
            labelstest = alllabels[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]
            
        if 5==6: # mix tel aviv and michaelstadt
            talllabelsu=np.load('tgridux.npy').ravel()
            talllabelsv=np.load('tgriduy.npy').ravel()
            talllabelsw=np.load('tgriduz.npy').ravel()
            talllabels = (talllabelsu**2+talllabelsv**2+talllabelsw**2)**0.5
            del talllabelsu, talllabelsv, talllabelsw
            
            tallfeatures=np.load('taf.npy')
            tgrid_y=np.load('tgrid_y.npy').ravel()
            talllabels=talllabels[::jpjp]
            tallfeatures=tallfeatures[::jpjp]
            tgrid_y=tgrid_y[::jpjp]
            
            # allfeatures=np.asarray(allfeatures)

            # tmsklabel = talllabels == 0.0
            # tfeatures0 = tallfeatures[~tmsklabel]
            # tlabels0 = talllabels[~tmsklabel]

            # msklabel = alllabels == 0.0
            # features0 = allfeatures[~msklabel]
            # labels0 = alllabels[~msklabel]
            
            scaler = preprocessing.MinMaxScaler()
            print('before scaler fit', datetime.datetime.now())
            scaler.fit(allfeatures)
            print('after scaler fit', datetime.datetime.now())
            strainx = scaler.transform(allfeatures)
            print('strainx', datetime.datetime.now())
            stestx = scaler.transform(tallfeatures)

            show=''
            model=models[0][1]

            model.fit(strainx,alllabelsu)
            print('after fit', datetime.datetime.now())
            predictvaluesx = model.predict(stestx)
            print('after predict', datetime.datetime.now())
#                corr = np.corrcoef(testy[::jumps], predictvaluesx)[0, 1]
            mul=1.64
            corr = sklearn.metrics.r2_score(predictvaluesx,talllabels*mul)
            print( 'corr is ',mul, corr, np.corrcoef(predictvaluesx, talllabels*mul)[1, 0], len(talllabels))
            ln = np.linspace(0,7,8)
            plt.figure()
            plt.scatter(talllabels*mul,predictvaluesx,s=1)
            plt.plot(ln,ln,'-r')
            plt.xlabel('RANS')
            plt.ylabel('ML')          
        
            
            

        
        msklabel = labels == 0.0
        features0 = features[~msklabel]
        labels0 = labels[~msklabel]
#        labels0limit=labels0[~msklabel]

        if featurestest is None:
            print ('random sample')
            self.msk = np.random.rand(len(labels0)) < 0.9
            lens=len(features0[self.msk])
            jumps=int(lens//200000) # 1 
            print ('chosen sample0000, jumps', lens, jumps)
            trainx = features0[self.msk][::jumps]
            trainy = labels0[self.msk][::jumps]
            lens=len(features0[~self.msk])
            jumps=int(lens//100000) # 1 
            print ('chosen sample11111, jumps', lens, jumps)            
            testx = features0[~self.msk][::jumps]
            testy = labels0[~self.msk][::jumps]
        else:
            lens=len(features0)
            jumps=int(lens//200000) # 1 
            print ('chosen sample, jumps', lens, jumps)
            msktest = labels == 0.0
            features1 = features[~msktest]
            labels1 = labels[~msktest]
#            labels1u = labelsu[~msktest]
#            labels1v = labelsv[~msktest]
            trainx = features0[::jumps]
            trainy = labels0[::jumps]
            testx = features1[::jumps]
            testy = labels1[::jumps]
        
        bestscore = -9999.9
        bestmodel = None
        bestscaler = None

#        print('scaler', trainy)
        scaler = preprocessing.MinMaxScaler()
        print('before scaler fit', datetime.datetime.now())
        scaler.fit(trainx)
        print('after scaler fit', datetime.datetime.now())
        strainx = scaler.transform(trainx)
        print('strainx', datetime.datetime.now())
        stestx = scaler.transform(testx)

        if 5==6:
#            np.save('f0',features)
#            np.save('f1',featurestest)
#            np.save('l0',labels)
#            np.save('l1',labeltest)
            features=np.load('mf0.npy')
            featurestest=np.load('mf1.npy')
            labels=np.load('ml0.npy')
            labeltest=np.load('ml1.npy')
            
            show=''
            model=models[0][1]
            speedlimit=5
            jp=198
            sc=[]
            scn=[]
            for jp in range(198,199):
#                model.fit(strainx[::jp],trainy[::jp])
                jp2=int(len(trainy[trainy<speedlimit])*jp/1.)
                jp1=np.random.choice(len(trainy[trainy<speedlimit]),size=jp2, replace=False)
                model.fit(strainx[trainy<speedlimit],trainy[trainy<speedlimit])
                print('after fit', datetime.datetime.now())
                predictvaluesx = model.predict(stestx[testy<speedlimit])
                print('after predict', datetime.datetime.now())
#                corr = np.corrcoef(testy[::jumps], predictvaluesx)[0, 1]
                corr = sklearn.metrics.r2_score(testy[testy<speedlimit],predictvaluesx)
                sc.append(corr)
                scn.append(jp2)
                print(jp, 'corr is ',corr, len(trainx[jp1]))
            plt.figure()
            plt.plot(scn,sc)
            plt.xlabel('training size')
            plt.ylabel('$R^2$')
            plt.figure()
            aaa1=testy[testy<speedlimit][::jumps]
            aaa2=predictvaluesx
            plt.scatter(aaa1[aaa1!=aaa2],aaa2[aaa1!=aaa2],s=1)
            plt.scatter(testy[testy<speedlimit][::jumps], predictvaluesx,s=1)
            
            
            plt.figure()
            plt.scatter(testy[testy<speedlimit], predictvaluesx, s=1, c='b')
            plt.scatter(testy[::jumps], labels1u, s=1, c='r')
            plt.scatter(testy[::jumps], labels1v, s=1, c='g')
            plt.scatter(testy[::jumps], np.arctan2(labels1v,labels1u)*180/np.pi, s=1, c='k')
            plt.xlabel('RANS [m/s]')
            plt.ylabel('Machine learning [m/s]')
            plt.title(show)
            plt.show()  
            np.save('testydown.npy',testy)
            np.save('predictvaluesxdown.npy',predictvaluesx)
            a=np.load('testyup.npy')
            b=np.load('predictvaluesxup.npy')
            plt.scatter(a, b, s=1, c='b')
            
            aaa = (trainx[:,0]<110) & (trainx[:,1]<110) & (trainx[:,0]>4) & (trainx[:,1]>4) & (trainx[:,4]<8) & (trainx[:,4]>1)
            np.argwhere(aaa==True) 
            trainx[474268]
            labels[474268]
            
            streetwide=10
            for i in range(streetwide):
                street = [
#                            28., #grid_z.ravel(),  #1
                            i, #dplim, #4
                            streetwide - i, #dnlim, #4
                            100., # dllim, #4
                            100., #drlim, #4
                            28, #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            i, #grid_dp.ravel(), #4
                            i/8., #grid_dp.ravel()/grid_z.ravel(), #4
                            streetwide - i, #grid_dn.ravel(), #5
                            streetwide, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            100, #grid_dl.ravel(),  #6
                            100, #grid_dr.ravel(),  #7
                            200, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            30.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            0.6 #angle #12
                            ]
                print(i, model.predict(scaler.transform([np.asarray(street)])))
            
            streetwide  = 30
            parcelheight = 12
            streetlength = 44
            for i in range(streetwide):
                street = [
#                            parcelheight, #grid_z.ravel(),  #1
                            300., #dplim, #4
                            300., #dnlim, #4
                            i, # dllim, #4
                            streetwide - i, #drlim, #4
                            parcelheight, #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            streetlength, #grid_dp.ravel(), #4
                            streetlength/parcelheight, #grid_dp.ravel()/grid_z.ravel(), #4
                            streetlength, #grid_dn.ravel(), #5
                            streetlength+streetlength, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            i, #grid_dl.ravel(),  #6
                            streetwide-i, #grid_dr.ravel(),  #7
                            streetwide, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            40.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            0.4 #angle #12
                            ]
#                print(i,street)
                
                print(i, model.predict(scaler.transform([np.asarray(street)])))

############ m feature 474268
                
            print(labels[474268])
            print(trainx[474268])
            streetwide=int(1.6949152e+01+4.2372880e+00)
#            streetwide=int(4.8694424e+01+4.2342978e+00)
            
            for i in range(streetwide):
                street = [
#                            28., #grid_z.ravel(),  #1
                            4.8694424e+01, #dplim, #4
                            4.2342978e+00, #dnlim, #4
                            100., # dllim, #4
                            100., #drlim, #4
                            28, #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            i, #grid_dp.ravel(), #4
                            i/8., #grid_dp.ravel()/grid_z.ravel(), #4
                            streetwide - i, #grid_dn.ravel(), #5
                            streetwide, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            100, #grid_dl.ravel(),  #6
                            100, #grid_dr.ravel(),  #7
                            200, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            30.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            0.6 #angle #12
                            ]
                
                street = np.asarray([4.8694424e+01, 4.2342978e+00, i, streetwide-i,
       2.1428571e+00, 1.0000000e+00, 4.8694424e+01, 2.2724066e+01,
       4.2342978e+00, 5.2928722e+01, i, streetwide-i,
       streetwide, 1.3973182e+03, 3.8074985e-01])

                street = np.asarray([4.8694424e+01, 4.2342978e+00, 1.6949152e+01, 4.2372880e+00,
       2.1428571e+00, 1.0000000e+00, 4.8694424e+01, 2.2724066e+01,
       4.2342978e+00, 5.2928722e+01, 1.6949152e+01, 4.2372880e+00,
       2.1186440e+01, 1.3973182e+03, 3.8074985e-01])
                
                print(i, model.predict(scaler.transform([np.asarray(street)])))
                
                


############  up the height
            streetwide  = 10
            streetlength = 220
            mlhigh=[]
            prhigh=[]
            for i in range(1,80):
                street = [
#                            parcelheight, #grid_z.ravel(),  #1
                            streetlength, #dplim, #4
                            streetlength, #dnlim, #4
                            streetwide, # dllim, #4
                            streetwide, #drlim, #4
                            i/2., #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            streetlength, #grid_dp.ravel(), #4
                            streetlength/(i/2.), #grid_dp.ravel()/grid_z.ravel(), #4
                            streetlength, #grid_dn.ravel(), #5
                            streetlength+streetlength, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            streetwide, #grid_dl.ravel(),  #6
                            streetwide, #grid_dr.ravel(),  #7
                            streetwide+streetwide, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            20.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            1.4 #angle #12
                            ]
                mlhigh.append(model.predict(scaler.transform([np.asarray(street)]))[0])
                prhigh.append(windprofile(i/2, uref=3, href=24, he=10, lambdap=0.3, lambdaf=0.8,beta=0.3, verbose=False))
#                print(i, model.predict(scaler.transform([np.asarray(street)])))

            plt.figure()
            zhigh=np.linspace(1,len(mlhigh), len(mlhigh))/2.
            plt.plot(mlhigh, zhigh,'b', label='machine learning')
            plt.plot(.24*np.log(zhigh/0.00002), zhigh,'r', label='logarithmic profile ')
#            plt.plot(prhigh, zhigh,'g', label='exp - log profile ')
            plt.xlabel('u [m/s]')
            plt.ylabel('z [m]')
            plt.legend()

        
        for name,model in self.models:
            if name !='KNN6' or show!='':            
                print('model:', name, datetime.datetime.now())
                model.fit(strainx,trainy)
                print('submit:', name, datetime.datetime.now())           
    #            print('after scale transform time', datetime.datetime.now(), type(stestx), stestx)
                predictvaluesx = model.predict(stestx)
                corr = np.corrcoef(testy, predictvaluesx)[0, 1]
                print('<<<ML>>>',name, corr, datetime.datetime.now())
                if corr>bestscore:
                        bestscore = corr
                        bestmodel = model
                        bestscaler = scaler

#        for name,model in self.models:
#            print('model:', name, datetime.datetime.now())
#            p = None
##            p = multiprocessing.Process(target=model.fit,args=(strainx,trainy))
#            p = multiprocessing.Process(target=model.fit,args=(trainx,trainy))
#            print('submit:', name, datetime.datetime.now())           
#            p.start()
#            print('started:', name, datetime.datetime.now())
#            p.join(20)
#            print('waited:', name, datetime.datetime.now())
#            if p.is_alive():
#                print ("running too long... let's kill it...", datetime.datetime.now())
#                p.terminate()
#                p.join()
#            else:
#                print('after reasonable time', datetime.datetime.now())
#                stestx = scaler.transform(testx)
##                print('mid',stestx[0:3])
#                print('after scale transform time', datetime.datetime.now(), type(stestx), stestx)
##                predictvaluesx = model.predict(stestx.T)
#                predictvaluesx = model.predict(testx)
#                print('before corr', datetime.datetime.now())  
#                corr = np.corrcoef(testy, predictvaluesx)[0, 1]
#                print('<<<ML>>>',name, corr, datetime.datetime.now())
#                if corr>bestscore:
#                    bestscore = corr
#                    bestmodel = model
#                    bestscaler = scaler

        self.model = bestmodel
        self.scaler = bestscaler
        self.score = bestscore
        print('test2')
        if show!='':
            testmin = testy.min()
            testmax = testy.max()
            test11 = np.linspace(testmin, testmax, 100)
            predictvaluesx = self.model.predict(stestx)
#            predictvaluesx = self.predict(testx)
            plt.figure()
            plt.scatter(testy, predictvaluesx, s=1, c='b')
            plt.scatter(test11, test11, s=1, c='r')
            plt.xlabel('testy')
            plt.ylabel('predictvalues')
            plt.title(show)
            plt.show()
#            plt.figure()
#            testxnp = np.asarray(testx)
#            sc = plt.scatter(testxnp[:, 0], testxnp[:, 1], c=(testy - predictvaluesx), s=3)  # , cmap=cm)
#            plt.colorbar(sc)
#            plt.xlabel('testxnp - X')
#            plt.ylabel('testxnp - Y')
#            plt.title(show)
#            plt.show()

#        predictfeatures = features[12345]            
#        test2 = self.predict(predictfeatures)
#        print ('1<<',features[12345])
#        print ('1>>',test2)
#
#        fea =  np.asarray([3.92786880e+00,    3.92786880e+00,   2.00e+00,   2.00000000e+00, 
#           4.00,   800.0e+00,   600.e+00,           8.52369517e-02])
#        modfea = self.predict(fea)
#        print ('>>111',modfea)        

        print ('fin fit', show, bestmodel, bestscore, datetime.datetime.now())
        
        return bestmodel, bestscaler, bestscore

    def info(self):
        print('model:', self.model)
        print('scaler:', self.scaler)
        print('score:', self.score)
        # return self.score()

    def predict(self, features):
        scaled = self.scaler.transform(features)
#        labelspredict = self.model.predict(scaler.transform(features))
#        labelspredict = self.model.predict(features)
        labelspredict = self.model.predict(scaled)
        return labelspredict


    def corr_vector(testyux, testyuy, testyuz, predictvaluesx, predictvaluesy, predictvaluesz):
        ttt = []
        fff = []
        ttt.append(np.asarray(testyux).ravel())
        ttt.append(np.asarray(testyuy).ravel())
        ttt.append(np.asarray(testyuz).ravel())
        nc = len(ttt)
        print('n_components=', nc)
        ttt = np.asarray(ttt).T
        print(ttt.shape)
        fff.append(np.asarray(predictvaluesx).ravel())
        fff.append(np.asarray(predictvaluesy).ravel())
        fff.append(np.asarray(predictvaluesz).ravel())
        fff = np.asarray(fff).T
        cca = CCA(n_components=nc)
        cca.fit(ttt, fff)
        result = cca.score(ttt, fff)
        print('cca score1', result)

        return result
     
        
#https://stackoverflow.com/questions/492519/timeout-on-a-function-call
# bar
#def bar():
#    for i in range(100):
#        print "Tick"
#        time.sleep(1)
#
#if __name__ == '__main__':
#    # Start bar as a process
#    p = multiprocessing.Process(target=bar)
#    p.start()
#
#    # Wait for 10 seconds or until process finishes
#    p.join(10)
#
#    # If thread is still active
#    if p.is_alive():
#        print "running... let's kill it..."
#
#        # Terminate
#        p.terminate()
#        p.join()

if __name__ == '__main__':
    print('main ml')
    gridztlv = np.load('/ibdata2/nirb/Projects/tlvz.npy')
    griduxtlv = np.load('/ibdata2/nirb/Projects/tlvux.npy')
    
    
    ml2 = ml()
    clf2, scaler, score = ml2.fit(features0, labelsu20, show='U3', featurestest = features1, labeltest = labelsu21)
    ml2.save('ml3-'+learnfile)

    u2ml = ml2.predict(features1)
